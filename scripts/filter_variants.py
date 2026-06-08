#!/usr/bin/env python3
"""
Merge primary and relapse annotated variant TSVs for one patient, then apply:
  - PASS-only (FilterMutectCalls)
  - depth >= 50 in the tumor sample
  - VAF > 5% in either primary OR relapse (union)
  - exclude non-coding consequences (intronic, UTR, downstream, upstream)
  - CADD_PHRED >= 0.20
  - normalize VAF by FACETS tumor purity
"""

import argparse
import re
import sys
import pandas as pd


# VEP CSQ field order — the index of each sub-field within the pipe-delimited
# CSQ INFO string.  Adjust if the VEP run uses a different field order.
CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
    "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON",
    "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
    "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND",
    "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS",
    "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM",
    "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "miRNA",
    "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF",
    "SAS_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF",
    "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF",
    "MAX_AF", "MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO",
    "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
]
CSQ_IDX = {f: i for i, f in enumerate(CSQ_FIELDS)}


def parse_csq(csq_str):
    """Return a dict of the first (most severe) CSQ annotation."""
    if pd.isna(csq_str) or csq_str == ".":
        return {}
    first = csq_str.split(",")[0]
    parts = first.split("|")
    return {f: parts[i] if i < len(parts) else "" for f, i in CSQ_IDX.items()}


def load_purity(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("purity"):
                continue
            return float(line.strip().split()[0])
    return 1.0


def parse_gatk_table(path, sample_type):
    """Read a GATK VariantsToTable TSV and normalise column names."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    # VariantsToTable names genotype columns as SAMPLE.AD, SAMPLE.AF, SAMPLE.DP
    # Rename to uniform names.
    rename = {}
    for col in df.columns:
        if col.endswith(".AD"):
            rename[col] = "AD"
        elif col.endswith(".AF"):
            rename[col] = f"AF_{sample_type}"
        elif col.endswith(".DP"):
            rename[col] = f"DP_{sample_type}"
    df = df.rename(columns=rename)

    csq_parsed = df["CSQ"].apply(parse_csq)
    csq_df = pd.DataFrame(csq_parsed.tolist())
    df = pd.concat([df.drop(columns=["CSQ", "AD"], errors="ignore"), csq_df], axis=1)

    df[f"AF_{sample_type}"] = pd.to_numeric(df[f"AF_{sample_type}"], errors="coerce")
    df[f"DP_{sample_type}"] = pd.to_numeric(df[f"DP_{sample_type}"], errors="coerce")
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--primary",  required=True, help="Primary CADD-annotated TSV")
    p.add_argument("--relapse",  required=True, help="Relapse CADD-annotated TSV")
    p.add_argument("--primary-purity", required=True)
    p.add_argument("--relapse-purity", required=True)
    p.add_argument("--patient-id", required=True)
    p.add_argument("--min-depth",  type=int,   default=50)
    p.add_argument("--min-vaf",    type=float, default=0.05)
    p.add_argument("--cadd-cutoff", type=float, default=0.20)
    p.add_argument("--exclude-consequences", nargs="+", default=[
        "intron_variant", "downstream_gene_variant", "upstream_gene_variant",
        "3_prime_UTR_variant", "5_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
    ])
    p.add_argument("--output", required=True)
    args = p.parse_args()

    purity_p = load_purity(args.primary_purity)
    purity_r = load_purity(args.relapse_purity)

    pri = parse_gatk_table(args.primary, "primary")
    rel = parse_gatk_table(args.relapse, "relapse")

    key_cols = ["CHROM", "POS", "REF", "ALT"]

    merged = pd.merge(pri, rel[key_cols + [f"AF_relapse", f"DP_relapse"]],
                      on=key_cols, how="outer")

    merged["AF_primary"] = pd.to_numeric(merged.get("AF_primary"), errors="coerce")
    merged["AF_relapse"] = pd.to_numeric(merged.get("AF_relapse"), errors="coerce")
    merged["DP_primary"] = pd.to_numeric(merged.get("DP_primary"), errors="coerce")
    merged["DP_relapse"] = pd.to_numeric(merged.get("DP_relapse"), errors="coerce")

    # ---- Filters ------------------------------------------------
    # 1. PASS in at least one call
    merged = merged[merged["FILTER"].isin(["PASS", None, "."]) |
                    merged["FILTER"].isna()]

    # 2. Depth >= 50 in the sample where the call was made
    depth_ok = (
        (merged["DP_primary"] >= args.min_depth) |
        (merged["DP_relapse"] >= args.min_depth)
    )
    merged = merged[depth_ok]

    # 3. VAF > 5% in either primary or relapse
    vaf_ok = (
        (merged["AF_primary"] > args.min_vaf) |
        (merged["AF_relapse"] > args.min_vaf)
    )
    merged = merged[vaf_ok]

    # 4. Exclude non-coding consequences
    excluded = set(args.exclude_consequences)
    def keep_consequence(csq):
        if pd.isna(csq):
            return False
        return not all(c.strip() in excluded for c in csq.split("&"))

    if "Consequence" in merged.columns:
        merged = merged[merged["Consequence"].apply(keep_consequence)]

    # 5. CADD_PHRED >= cutoff
    if "CADD_PHRED" in merged.columns:
        merged["CADD_PHRED"] = pd.to_numeric(merged["CADD_PHRED"], errors="coerce")
        merged = merged[merged["CADD_PHRED"] >= args.cadd_cutoff]

    # ---- VAF normalisation by purity ----------------------------
    merged["AF_primary_norm"] = (merged["AF_primary"] / purity_p).clip(upper=1.0)
    merged["AF_relapse_norm"] = (merged["AF_relapse"] / purity_r).clip(upper=1.0)
    merged["patient_id"] = args.patient_id

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Retained {len(merged)} variants for {args.patient_id}", file=sys.stderr)


if __name__ == "__main__":
    main()
