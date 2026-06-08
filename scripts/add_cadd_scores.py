#!/usr/bin/env python3
"""
Look up CADD v1.6 pre-scored files (tabix-indexed) for each variant in a
GATK VariantsToTable TSV and append CADD_RAW and CADD_PHRED columns.

SNV scores come from whole_genome_SNVs.tsv.gz.
InDel scores come from InDels.tsv.gz.

Both files have the format (1-based, tab-separated, bgzip + tabix):
  #Chrom  Pos  Ref  Alt  RawScore  PHRED
"""

import argparse
import pysam
import pandas as pd
import sys


def strip_chr(chrom):
    """Normalise chromosome name to match CADD file (no 'chr' prefix)."""
    return chrom.replace("chr", "").replace("Chr", "")


def lookup_cadd(tabix_file, chrom, pos, ref, alt):
    """Return (raw, phred) tuple from a tabix-indexed CADD file, or (None, None)."""
    chrom = strip_chr(chrom)
    try:
        pos = int(pos)
        for row in tabix_file.fetch(chrom, pos - 1, pos):
            fields = row.split("\t")
            if fields[2] == ref and fields[3] == alt:
                return float(fields[4]), float(fields[5])
    except (ValueError, KeyError, pysam.TabixError):
        pass
    return None, None


def is_indel(ref, alt):
    return len(ref) != len(alt)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",      required=True, help="GATK VariantsToTable TSV")
    p.add_argument("--snv-file",   required=True, help="CADD SNV .tsv.gz (tabix-indexed)")
    p.add_argument("--indel-file", required=True, help="CADD InDel .tsv.gz (tabix-indexed)")
    p.add_argument("--output",     required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    snv_tbx   = pysam.TabixFile(args.snv_file)
    indel_tbx = pysam.TabixFile(args.indel_file)

    raw_scores  = []
    phred_scores = []

    for _, row in df.iterrows():
        chrom = row.get("CHROM", row.get("#CHROM", ""))
        pos   = row.get("POS", "")
        ref   = row.get("REF", "")
        alt   = row.get("ALT", "")

        tbx = indel_tbx if is_indel(ref, alt) else snv_tbx
        raw, phred = lookup_cadd(tbx, chrom, pos, ref, alt)
        raw_scores.append(raw)
        phred_scores.append(phred)

    snv_tbx.close()
    indel_tbx.close()

    df["CADD_RAW"]   = raw_scores
    df["CADD_PHRED"] = phred_scores

    missing = df["CADD_PHRED"].isna().sum()
    if missing:
        print(f"WARNING: {missing}/{len(df)} variants had no CADD score.", file=sys.stderr)

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
