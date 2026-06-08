#!/usr/bin/env python3
"""
Mutual exclusivity analysis using the DISCOVER Python package.
https://github.com/NKI-CCB/DISCOVER

The event matrix treats primary and relapse as separate columns
(matching the approach in 1554.py: pivot on patientID × SampleType).
Only genes present in > min_samples columns are tested.

Install: conda install -c https://ccb.nki.nl/software/discover/repos/conda discover

Output TSV columns:
  gene_a, gene_b, n_a, n_b, n_both, p_value, significant (p < 0.05)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

import discover


def build_event_matrix(results_dir, patients, min_samples):
    """
    Build a binary gene × (patient, sample_type) event matrix.

    Columns are a MultiIndex of (patient_id, sample_type) — Primary and Relapse
    are kept as separate columns, consistent with the project's existing analysis.
    A cell is 1 if AF > 0 in that sample after all filters have been applied.

    Returns a pandas DataFrame with MultiIndex columns, dtype float (0/1).
    """
    records = []  # list of (patient_id, sample_type, gene) tuples

    for pid in patients:
        path = os.path.join(results_dir, pid, "final_variants", f"{pid}_final.tsv")
        if not os.path.exists(path):
            print(f"WARNING: missing {path}", file=sys.stderr)
            continue

        df = pd.read_csv(path, sep="\t")
        gene_col = next((c for c in ("SYMBOL", "Gene", "gene") if c in df.columns), None)
        if gene_col is None:
            print(f"WARNING: no gene column in {path}", file=sys.stderr)
            continue

        df["AF_primary"] = pd.to_numeric(df.get("AF_primary"), errors="coerce").fillna(0)
        df["AF_relapse"] = pd.to_numeric(df.get("AF_relapse"), errors="coerce").fillna(0)

        for gene in df.loc[df["AF_primary"] > 0, gene_col].dropna().unique():
            records.append((pid, "Primary", gene))
        for gene in df.loc[df["AF_relapse"] > 0, gene_col].dropna().unique():
            records.append((pid, "Relapse", gene))

    if not records:
        raise ValueError("No mutation records found. Check results_dir and patient IDs.")

    long = pd.DataFrame(records, columns=["patientID", "SampleType", "SYMBOL"])
    long["presence"] = 1

    # Deduplicate (a gene can appear multiple times per sample from different mutations)
    long = long.drop_duplicates(subset=["patientID", "SampleType", "SYMBOL"])

    mat = long.pivot(
        index="SYMBOL",
        columns=["patientID", "SampleType"],
        values="presence",
    ).fillna(0)

    # Keep genes present in > min_samples columns (same filter as 1554.py: sum > 4)
    keep = mat.sum(axis=1) > min_samples
    mat = mat.loc[keep]

    if len(mat) < 2:
        raise ValueError(
            f"Only {len(mat)} gene(s) present in > {min_samples} sample columns. "
            "Lower --min-samples or check input files."
        )

    return mat


def flatten_pvalues(pval_matrix, mat):
    """
    Convert the upper triangle of a gene × gene p-value matrix to a long DataFrame.
    Adds observed counts (n_a, n_b, n_both) from the event matrix.
    """
    genes = pval_matrix.index.tolist()
    rows = []
    for i, ga in enumerate(genes):
        for j, gb in enumerate(genes):
            if j <= i:
                continue
            pval = pval_matrix.iloc[i, j]
            if pd.isna(pval):
                continue
            n_a    = int(mat.loc[ga].sum())
            n_b    = int(mat.loc[gb].sum())
            n_both = int((mat.loc[ga].values.astype(bool) &
                          mat.loc[gb].values.astype(bool)).sum())
            rows.append({
                "gene_a":      ga,
                "gene_b":      gb,
                "n_a":         n_a,
                "n_b":         n_b,
                "n_both":      n_both,
                "p_value":     float(pval),
                "significant": pval < 0.05,
            })

    return pd.DataFrame(rows).sort_values("p_value")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--patients",    nargs="+", required=True)
    p.add_argument("--min-samples", type=int, default=4,
                   help="Genes in > N sample-columns are tested (default 4, matching sum>4 filter)")
    p.add_argument("--output",      required=True)
    args = p.parse_args()

    # ---- Build binary event matrix -----------------------------------------
    mat = build_event_matrix(args.results_dir, args.patients, args.min_samples)
    print(
        f"Event matrix: {mat.shape[0]} genes × {mat.shape[1]} sample-columns "
        f"({len(args.patients)} patients × primary/relapse)",
        file=sys.stderr,
    )

    # ---- Fit DISCOVER background model and run mutual exclusivity test -----
    events = discover.DiscoverMatrix(mat)
    result_mutex = discover.pairwise_discover_test(events)   # default: alternative="less"

    # ---- Flatten p-value matrix to long TSV --------------------------------
    out = flatten_pvalues(result_mutex.pvalues, mat)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    n_sig = out["significant"].sum()
    print(
        f"DISCOVER: {len(out)} gene pairs tested; {n_sig} significant (p < 0.05)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
