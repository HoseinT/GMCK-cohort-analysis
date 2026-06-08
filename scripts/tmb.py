#!/usr/bin/env python3
"""
Tumor Mutational Burden (TMB) per sample.

For each patient, the final filtered variant TSV already has one row per
unique variant position, with AF_primary and AF_relapse columns. We count
variants that are present in the primary sample and variants present in the
relapse sample separately, then divide by the panel size in megabases.

A variant is counted for a sample type if its AF in that sample is > 0
(i.e., it was observed in that sample after filtering).
"""

import argparse
import os
import sys
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir",   required=True)
    p.add_argument("--patients",      nargs="+", required=True)
    p.add_argument("--panel-size-mb", type=float, default=2.4)
    p.add_argument("--output",        required=True)
    args = p.parse_args()

    records = []

    for patient in args.patients:
        tsv_path = os.path.join(
            args.results_dir, patient, "final_variants", f"{patient}_final.tsv"
        )
        if not os.path.exists(tsv_path):
            print(f"WARNING: missing {tsv_path}", file=sys.stderr)
            continue

        df = pd.read_csv(tsv_path, sep="\t")
        df["AF_primary"] = pd.to_numeric(df.get("AF_primary"), errors="coerce").fillna(0)
        df["AF_relapse"] = pd.to_numeric(df.get("AF_relapse"), errors="coerce").fillna(0)

        for sample_type, af_col in [("primary", "AF_primary"), ("relapse", "AF_relapse")]:
            n_snvs = int((df[af_col] > 0).sum())
            tmb    = n_snvs / args.panel_size_mb
            records.append({
                "patient_id":   patient,
                "sample_type":  sample_type,
                "n_snvs":       n_snvs,
                "tmb":          round(tmb, 4),
                "panel_size_mb": args.panel_size_mb,
            })

    result = pd.DataFrame(records)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"TMB written for {len(args.patients)} patients.", file=sys.stderr)


if __name__ == "__main__":
    main()
