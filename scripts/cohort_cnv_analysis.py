#!/usr/bin/env python3
"""
Cohort-level ploidy-adjusted CNV analysis.

Steps:
  1. Map FACETS segments to panel genes (BED4 intersection).
  2. Compute adjusted TCN = tcn.em - ploidy for every gene × sample.
  3. Per-sample CNV burden: count amp (adj_tcn > threshold) and del (adj_tcn < threshold).
  4. Recurrent genes: mean adj_tcn across samples vs avg thresholds.
  5. Paired Wilcoxon signed-rank test (primary vs relapse) with BH correction.
  6. Write adjusted TCN matrix, burden, recurrent, and Wilcoxon result tables.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


# ---- Segment → gene mapping -------------------------------------------------

def load_gene_bed(bed_path):
    """Load BED4 gene coordinates into a DataFrame."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "start", "end", "gene"],
        dtype={"chrom": str, "start": int, "end": int, "gene": str},
    )
    df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)
    return df


def assign_gene_tcn(segments, gene_bed):
    """
    For each gene, assign the adjusted TCN from the overlapping FACETS segment
    with the largest overlap (in bp). Returns a dict {gene: adj_tcn}.
    """
    seg = segments.copy()
    seg["chrom"] = seg["chrom"].astype(str).str.replace("chr", "", regex=False)

    result = {}
    for _, gene_row in gene_bed.iterrows():
        gene  = gene_row["gene"]
        chrom = gene_row["chrom"]
        gstart = gene_row["start"]
        gend   = gene_row["end"]

        candidates = seg[
            (seg["chrom"] == chrom) &
            (seg["start"] <= gend) &
            (seg["end"]   >= gstart)
        ].copy()

        if candidates.empty:
            result[gene] = np.nan
            continue

        # Overlap length as tiebreaker
        candidates["overlap"] = (
            candidates["end"].clip(upper=gend) -
            candidates["start"].clip(lower=gstart)
        )
        best = candidates.loc[candidates["overlap"].idxmax()]
        result[gene] = best["adj_tcn"]

    return result


def load_purity(purity_path):
    """Return (purity, ploidy) from a FACETS purity file."""
    with open(purity_path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    # Header line: "purity\tploidy", data line: "0.85\t2.1"
    for line in lines:
        parts = line.split()
        if parts[0] in ("purity", "NA"):
            continue
        try:
            return float(parts[0]), float(parts[1])
        except (ValueError, IndexError):
            continue
    return np.nan, np.nan


# ---- Main -------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir",           required=True)
    p.add_argument("--patients",              nargs="+", required=True)
    p.add_argument("--gene-bed",              required=True)
    p.add_argument("--amp-threshold",         type=float, default=1.0)
    p.add_argument("--del-threshold",         type=float, default=-1.0)
    p.add_argument("--avg-amp-threshold",     type=float, default=0.85)
    p.add_argument("--avg-del-threshold",     type=float, default=-0.50)
    p.add_argument("--burden-abs",            type=float, default=1.0)
    p.add_argument("--wilcoxon-min-patients", type=int,   default=5)
    p.add_argument("--outdir",                required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gene_bed = load_gene_bed(args.gene_bed)

    # ---- Build adjusted TCN matrix -----------------------------------------
    # Columns: {patient}_primary, {patient}_relapse
    adj_tcn_dict = {}

    for patient in args.patients:
        for tumor in ("primary", "relapse"):
            cnv_path    = os.path.join(args.results_dir, patient, "facets", f"{tumor}_cnv.tsv")
            purity_path = os.path.join(args.results_dir, patient, "facets", f"{tumor}_purity.txt")

            if not os.path.exists(cnv_path) or not os.path.exists(purity_path):
                print(f"WARNING: missing FACETS output for {patient}/{tumor}", file=sys.stderr)
                continue

            _, ploidy = load_purity(purity_path)
            seg = pd.read_csv(cnv_path, sep="\t")

            # Normalise column names (FACETS may output chrom/start/end or loc.start/loc.end)
            seg = seg.rename(columns={
                "loc.start": "start",
                "loc.end":   "end",
                "tcn.em":    "tcn_em",
            })
            if "chrom" not in seg.columns and "Chromosome" in seg.columns:
                seg = seg.rename(columns={"Chromosome": "chrom"})
            if "start" not in seg.columns and "Start" in seg.columns:
                seg = seg.rename(columns={"Start": "start", "End": "end"})

            seg["adj_tcn"] = seg["tcn_em"].astype(float) - float(ploidy)

            col_name = f"{patient}_{tumor}"
            adj_tcn_dict[col_name] = assign_gene_tcn(seg, gene_bed)

    matrix = pd.DataFrame(adj_tcn_dict, index=gene_bed["gene"].unique())
    matrix.index.name = "gene"
    matrix.to_csv(os.path.join(args.outdir, "cnv_adjusted_matrix.tsv"), sep="\t")

    # ---- Per-sample CNV burden ---------------------------------------------
    burden_records = []
    for col in matrix.columns:
        patient, tumor = col.rsplit("_", 1)
        vals = matrix[col].dropna()
        n_amp = int((vals > args.amp_threshold).sum())
        n_del = int((vals < args.del_threshold).sum())
        burden_records.append({
            "patient_id":    patient,
            "sample_type":   tumor,
            "n_amp":         n_amp,
            "n_del":         n_del,
            "total_burden":  n_amp + n_del,
        })

    burden_df = pd.DataFrame(burden_records)

    # Median and IQR per sample type
    for stype in ("primary", "relapse"):
        sub = burden_df[burden_df["sample_type"] == stype]["total_burden"]
        if not sub.empty:
            print(
                f"{stype} burden — median: {sub.median():.1f}, "
                f"Q1: {sub.quantile(0.25):.1f}, Q3: {sub.quantile(0.75):.1f}",
                file=sys.stderr,
            )

    burden_df.to_csv(os.path.join(args.outdir, "cnv_burden.tsv"), sep="\t", index=False)

    # ---- Recurrent CNV genes -----------------------------------------------
    primary_cols = [c for c in matrix.columns if c.endswith("_primary")]
    relapse_cols = [c for c in matrix.columns if c.endswith("_relapse")]

    recurrent = pd.DataFrame(index=matrix.index)
    recurrent["mean_adj_tcn_primary"] = matrix[primary_cols].mean(axis=1)
    recurrent["mean_adj_tcn_relapse"] = matrix[relapse_cols].mean(axis=1)
    recurrent["n_samples_primary"] = len(primary_cols)
    recurrent["n_samples_relapse"] = len(relapse_cols)
    recurrent["n_amp_primary"] = (matrix[primary_cols] > args.amp_threshold).sum(axis=1)
    recurrent["n_del_primary"] = (matrix[primary_cols] < args.del_threshold).sum(axis=1)
    recurrent["n_amp_relapse"] = (matrix[relapse_cols] > args.amp_threshold).sum(axis=1)
    recurrent["n_del_relapse"] = (matrix[relapse_cols] < args.del_threshold).sum(axis=1)

    recurrent["recurrent_amp_primary"] = recurrent["mean_adj_tcn_primary"] >= args.avg_amp_threshold
    recurrent["recurrent_del_primary"] = recurrent["mean_adj_tcn_primary"] <= args.avg_del_threshold
    recurrent["recurrent_amp_relapse"] = recurrent["mean_adj_tcn_relapse"] >= args.avg_amp_threshold
    recurrent["recurrent_del_relapse"] = recurrent["mean_adj_tcn_relapse"] <= args.avg_del_threshold

    recurrent.to_csv(os.path.join(args.outdir, "cnv_recurrent.tsv"), sep="\t")

    # ---- Paired Wilcoxon signed-rank test ----------------------------------
    # Match patients that have both primary and relapse columns.
    shared_patients = [
        pid for pid in args.patients
        if f"{pid}_primary" in matrix.columns and f"{pid}_relapse" in matrix.columns
    ]

    wil_records = []
    for gene in matrix.index:
        pri_vals = matrix.loc[gene, [f"{p}_primary" for p in shared_patients]]
        rel_vals = matrix.loc[gene, [f"{p}_relapse" for p in shared_patients]]

        # Drop pairs where either value is NaN
        mask = pri_vals.notna() & rel_vals.notna()
        pri_vals = pri_vals[mask].values.astype(float)
        rel_vals = rel_vals[mask].values.astype(float)

        if len(pri_vals) < args.wilcoxon_min_patients:
            continue

        # Wilcoxon requires non-zero differences; handle zero-difference gracefully
        if np.all(pri_vals == rel_vals):
            wil_records.append({
                "gene": gene,
                "mean_primary": float(np.mean(pri_vals)),
                "mean_relapse": float(np.mean(rel_vals)),
                "statistic": np.nan,
                "p_value": 1.0,
                "n_pairs": int(len(pri_vals)),
            })
            continue

        try:
            stat, pval = wilcoxon(pri_vals, rel_vals, zero_method="wilcox")
        except ValueError:
            stat, pval = np.nan, np.nan

        wil_records.append({
            "gene":         gene,
            "mean_primary": float(np.mean(pri_vals)),
            "mean_relapse": float(np.mean(rel_vals)),
            "statistic":    float(stat) if stat is not np.nan else np.nan,
            "p_value":      float(pval) if pval is not np.nan else np.nan,
            "n_pairs":      int(len(pri_vals)),
        })

    wil_df = pd.DataFrame(wil_records)
    if not wil_df.empty:
        valid = wil_df["p_value"].notna()
        fdr = np.full(len(wil_df), np.nan)
        if valid.sum() > 0:
            _, fdr_vals, _, _ = multipletests(
                wil_df.loc[valid, "p_value"].values, method="fdr_bh"
            )
            fdr[valid] = fdr_vals
        wil_df["fdr"] = fdr
        wil_df["delta_mean"] = wil_df["mean_relapse"] - wil_df["mean_primary"]
        wil_df = wil_df.sort_values("p_value")

    wil_df.to_csv(os.path.join(args.outdir, "cnv_wilcoxon.tsv"), sep="\t", index=False)
    print(f"Wilcoxon test run on {len(wil_df)} genes.", file=sys.stderr)


if __name__ == "__main__":
    main()
