#!/usr/bin/env python3
"""
Cohort-level Plotly visualisations.

Produces five HTML reports (+ optional PDF via kaleido):
  tmb.html                   — TMB grouped bar chart
  cnv_burden.html            — CNV burden stacked bar
  cnv_heatmap.html           — Ploidy-adjusted TCN heatmap (genes × samples)
  mutual_exclusivity.html    — DISCOVER gene-pair heatmap
  cnv_wilcoxon_volcano.html  — Volcano: Δmean adj-TCN vs -log10(FDR)
"""

import argparse
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def _save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.write_html(path)
    try:
        import kaleido  # noqa: F401
        fig.write_image(path.replace(".html", ".pdf"), engine="kaleido")
    except Exception:
        pass


# ---- TMB bar chart ---------------------------------------------------------

def plot_tmb(tmb_path, outdir):
    df = pd.read_csv(tmb_path, sep="\t")
    patients = df["patient_id"].unique().tolist()

    fig = go.Figure()
    colors = {"primary": "#2166ac", "relapse": "#d6604d"}

    for stype in ("primary", "relapse"):
        sub = df[df["sample_type"] == stype].set_index("patient_id")
        tmbs = [sub.loc[p, "tmb"] if p in sub.index else 0 for p in patients]
        fig.add_trace(go.Bar(
            name=stype.capitalize(),
            x=patients,
            y=tmbs,
            marker_color=colors[stype],
            text=[f"{v:.2f}" for v in tmbs],
            textposition="outside",
        ))

    fig.update_layout(
        title="Tumor Mutational Burden (SNVs / Mb)",
        xaxis_title="Patient",
        yaxis_title="TMB (mutations / Mb)",
        barmode="group",
        template="plotly_white",
        height=450,
        legend_title="Sample type",
    )
    _save(fig, os.path.join(outdir, "tmb.html"))


# ---- CNV burden stacked bar ------------------------------------------------

def plot_cnv_burden(burden_path, outdir):
    df = pd.read_csv(burden_path, sep="\t")
    df = df.sort_values("total_burden", ascending=False)
    df["label"] = df["patient_id"] + "_" + df["sample_type"]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="Amplifications",
        x=df["label"],
        y=df["n_amp"],
        marker_color="#a50f15",
    ))
    fig.add_trace(go.Bar(
        name="Deletions",
        x=df["label"],
        y=df["n_del"],
        marker_color="#1a1a8c",
    ))

    fig.update_layout(
        title="CNV Burden per Sample (|adjusted TCN| > 1)",
        xaxis_title="Sample",
        yaxis_title="Number of genes",
        barmode="stack",
        template="plotly_white",
        height=450,
        xaxis_tickangle=45,
    )
    _save(fig, os.path.join(outdir, "cnv_burden.html"))


# ---- Adjusted CNV heatmap --------------------------------------------------

def plot_cnv_heatmap(matrix_path, outdir):
    mat = pd.read_csv(matrix_path, sep="\t", index_col=0)

    # Keep genes with any non-neutral value in at least one sample
    mask = ((mat > 0.5) | (mat < -0.5)).any(axis=1)
    mat = mat[mask]
    if mat.empty:
        return

    fig = go.Figure(data=go.Heatmap(
        z=mat.values,
        x=mat.columns.tolist(),
        y=mat.index.tolist(),
        colorscale=[
            [0.0,  "#1a1a8c"],   # deep deletion
            [0.3,  "#6baed6"],   # mild deletion
            [0.5,  "#f7f7f7"],   # neutral
            [0.7,  "#fd8d3c"],   # gain
            [1.0,  "#a50f15"],   # amplification
        ],
        zmid=0,
        colorbar=dict(title="Adj. TCN"),
        hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Adj TCN: %{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="Ploidy-Adjusted Copy Number (genes × samples)",
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9)),
        template="plotly_white",
        height=max(400, 12 * len(mat)),
    )
    _save(fig, os.path.join(outdir, "cnv_heatmap.html"))


# ---- Mutual exclusivity heatmap --------------------------------------------

def plot_mutual_exclusivity(mutex_path, outdir):
    df = pd.read_csv(mutex_path, sep="\t")
    if df.empty:
        return

    # Only show genes involved in at least one significant pair
    sig = df[df["p_value"] < 0.05]
    if sig.empty:
        sig = df.nsmallest(min(50, len(df)), "p_value")

    genes = sorted(set(sig["gene_a"].tolist() + sig["gene_b"].tolist()))
    n = len(genes)
    if n < 2:
        return
    idx = {g: i for i, g in enumerate(genes)}

    mat = np.full((n, n), np.nan)
    hover = [["" for _ in genes] for _ in genes]

    for _, row in df[df["gene_a"].isin(genes) & df["gene_b"].isin(genes)].iterrows():
        i, j = idx[row["gene_a"]], idx[row["gene_b"]]
        val = -np.log10(max(row["p_value"], 1e-300))
        mat[i, j] = val
        mat[j, i] = val
        text = (
            f"{row['gene_a']} × {row['gene_b']}<br>"
            f"n_a={row['n_a']}  n_b={row['n_b']}  n_both={row['n_both']}<br>"
            f"p={row['p_value']:.3e}"
        )
        hover[i][j] = text
        hover[j][i] = text

    fig = go.Figure(data=go.Heatmap(
        z=mat,
        x=genes,
        y=genes,
        colorscale="Blues",
        colorbar=dict(title="−log10(p)"),
        text=hover,
        hoverinfo="text",
    ))
    fig.update_layout(
        title="Mutual Exclusivity — DISCOVER (significant pairs, p < 0.05)",
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9)),
        template="plotly_white",
        height=max(400, 18 * n),
        width=max(400, 18 * n),
    )
    _save(fig, os.path.join(outdir, "mutual_exclusivity.html"))


# ---- Wilcoxon volcano plot -------------------------------------------------

def plot_wilcoxon_volcano(wilcoxon_path, outdir):
    df = pd.read_csv(wilcoxon_path, sep="\t")
    df = df.dropna(subset=["fdr", "delta_mean"])

    df["neg_log_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))
    sig = df["fdr"] < 0.05

    colors = np.where(sig & (df["delta_mean"] > 0), "#a50f15",
             np.where(sig & (df["delta_mean"] < 0), "#1a1a8c", "#aaaaaa"))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["delta_mean"],
        y=df["neg_log_fdr"],
        mode="markers",
        marker=dict(color=colors.tolist(), size=7, opacity=0.8),
        text=df["gene"],
        hovertemplate="<b>%{text}</b><br>Δmean: %{x:.3f}<br>−log10(FDR): %{y:.2f}<extra></extra>",
        showlegend=False,
    ))

    # Label significant genes
    for _, row in df[sig].iterrows():
        fig.add_annotation(
            x=row["delta_mean"],
            y=row["neg_log_fdr"],
            text=row["gene"],
            showarrow=False,
            font=dict(size=8),
            xanchor="left",
        )

    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey",
                  annotation_text="FDR 0.05")
    fig.add_vline(x=0, line_dash="dot", line_color="grey")

    fig.update_layout(
        title="Paired Wilcoxon (primary vs relapse adj. TCN) — BH corrected",
        xaxis_title="Δ mean adj. TCN (relapse − primary)",
        yaxis_title="−log10(FDR)",
        template="plotly_white",
        height=500,
    )
    _save(fig, os.path.join(outdir, "cnv_wilcoxon_volcano.html"))


# ---- Entry point -----------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tmb",       required=True)
    p.add_argument("--burden",    required=True)
    p.add_argument("--recurrent", required=True)
    p.add_argument("--wilcoxon",  required=True)
    p.add_argument("--matrix",    required=True)
    p.add_argument("--mutex",     required=True)
    p.add_argument("--outdir",    required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    plot_tmb(args.tmb, args.outdir)
    plot_cnv_burden(args.burden, args.outdir)
    plot_cnv_heatmap(args.matrix, args.outdir)
    plot_mutual_exclusivity(args.mutex, args.outdir)
    plot_wilcoxon_volcano(args.wilcoxon, args.outdir)

    print("All cohort plots written to", args.outdir)


if __name__ == "__main__":
    main()
