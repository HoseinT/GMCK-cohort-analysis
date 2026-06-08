#!/usr/bin/env python3
"""
Two sub-commands:
  snv  — lollipop plot of primary vs relapse VAF for one patient
  cnv  — genome-wide CNV segment plot (primary + relapse) for one patient

Output is an interactive Plotly HTML file (and optionally a PDF via kaleido).

Usage:
  python visualize.py snv --variants <tsv> --patient-id <id> --output <html>
  python visualize.py cnv --primary-cnv <tsv> --relapse-cnv <tsv> \
                          --patient-id <id> --output <html>
"""

import argparse
import sys
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

CNV_COLORS = {
    "deep_loss":     "#1a1a8c",
    "loss":          "#6baed6",
    "CNLOH":         "#fdae6b",
    "neutral":       "#d3d3d3",
    "gain":          "#fd8d3c",
    "amplification": "#a50f15",
    "unknown":       "#cccccc",
}

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]


# ---- SNV lollipop --------------------------------------------------

def plot_snv(args):
    df = pd.read_csv(args.variants, sep="\t")
    if df.empty:
        fig = go.Figure()
        fig.update_layout(title=f"{args.patient_id} — no variants after filtering")
        fig.write_html(args.output)
        return

    df["label"] = df.apply(
        lambda r: f"{r.get('SYMBOL','?')} {r.get('HGVSp', r.get('HGVSc',''))}",
        axis=1,
    )
    df["AF_primary"] = pd.to_numeric(df.get("AF_primary"), errors="coerce").fillna(0)
    df["AF_relapse"] = pd.to_numeric(df.get("AF_relapse"), errors="coerce").fillna(0)
    df["AF_primary_norm"] = pd.to_numeric(df.get("AF_primary_norm"), errors="coerce").fillna(0)
    df["AF_relapse_norm"] = pd.to_numeric(df.get("AF_relapse_norm"), errors="coerce").fillna(0)

    hover = df.apply(
        lambda r: (
            f"<b>{r['label']}</b><br>"
            f"Position: {r.get('CHROM','')}:{r.get('POS','')}<br>"
            f"Change: {r.get('REF','')}>{r.get('ALT','')}<br>"
            f"Consequence: {r.get('Consequence','')}<br>"
            f"CADD PHRED: {r.get('CADD_PHRED','N/A')}<br>"
            f"Primary VAF: {r['AF_primary']:.3f} (norm: {r['AF_primary_norm']:.3f})<br>"
            f"Relapse VAF: {r['AF_relapse']:.3f} (norm: {r['AF_relapse_norm']:.3f})"
        ),
        axis=1,
    )

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        subplot_titles=["Primary VAF", "Relapse VAF"],
                        vertical_spacing=0.08)

    x = list(range(len(df)))

    for row_idx, (col, norm_col, title) in enumerate(
        [("AF_primary", "AF_primary_norm", "Primary"),
         ("AF_relapse", "AF_relapse_norm", "Relapse")],
        start=1,
    ):
        fig.add_trace(
            go.Scatter(
                x=x, y=df[col],
                mode="markers",
                marker=dict(size=10, color="#2166ac"),
                name=f"{title} VAF (raw)",
                text=hover,
                hoverinfo="text",
            ),
            row=row_idx, col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=x, y=df[norm_col],
                mode="markers",
                marker=dict(size=8, color="#d6604d", symbol="diamond"),
                name=f"{title} VAF (purity-norm)",
                text=hover,
                hoverinfo="text",
            ),
            row=row_idx, col=1,
        )
        for xi, yi in zip(x, df[col]):
            fig.add_shape(
                type="line",
                x0=xi, x1=xi, y0=0, y1=yi,
                line=dict(color="#2166ac", width=1),
                row=row_idx, col=1,
            )

    fig.update_xaxes(
        tickvals=x,
        ticktext=df["label"].tolist(),
        tickangle=45,
        row=2, col=1,
    )
    fig.update_yaxes(title_text="VAF", range=[0, 1.05])
    fig.update_layout(
        title=f"{args.patient_id} — Somatic Variants",
        height=700,
        showlegend=True,
        template="plotly_white",
    )

    fig.write_html(args.output)
    _try_pdf(fig, args.output)


# ---- CNV genome plot -----------------------------------------------

def _chr_to_int(chrom):
    c = str(chrom).replace("chr", "").replace("Chr", "")
    try:
        return int(c)
    except ValueError:
        return {"X": 23, "Y": 24}.get(c, 25)


def _build_offset_map(segs):
    """Compute cumulative base-pair offsets per chromosome for linear genome plot."""
    chr_max = {}
    for _, row in segs.iterrows():
        c = str(row["chrom"])
        end = int(row.get("end", row.get("loc.end", 0)))
        chr_max[c] = max(chr_max.get(c, 0), end)

    offsets = {}
    cumulative = 0
    for c in sorted(chr_max, key=_chr_to_int):
        offsets[c] = cumulative
        cumulative += chr_max[c]
    return offsets, cumulative


def plot_cnv(args):
    pri = pd.read_csv(args.primary_cnv, sep="\t")
    rel = pd.read_csv(args.relapse_cnv, sep="\t")

    for df in [pri, rel]:
        df["chrom"] = df["chrom"].astype(str).str.replace("chr", "")

    all_segs = pd.concat([pri, rel])
    offsets, total_len = _build_offset_map(all_segs)

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        subplot_titles=["Primary — Copy Number", "Relapse — Copy Number"],
        vertical_spacing=0.08,
    )

    def add_cnv_traces(df, row_idx):
        for _, seg in df.iterrows():
            chrom = str(seg["chrom"])
            start = int(seg.get("start", seg.get("loc.start", 0))) + offsets.get(chrom, 0)
            end   = int(seg.get("end",   seg.get("loc.end",   0))) + offsets.get(chrom, 0)
            cls   = seg.get("cnv_class", "unknown")
            tcn   = seg.get("tcn.em", seg.get("tcn", "?"))
            lcn   = seg.get("lcn.em", seg.get("lcn", "?"))
            color = CNV_COLORS.get(cls, "#cccccc")
            fig.add_shape(
                type="rect",
                x0=start, x1=end,
                y0=0, y1=1,
                fillcolor=color,
                line=dict(width=0),
                opacity=0.8,
                row=row_idx, col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=[(start + end) / 2],
                    y=[0.5],
                    mode="markers",
                    marker=dict(size=0.1, color=color),
                    hovertext=(
                        f"Chr{chrom}:{seg.get('start', seg.get('loc.start',''))}–"
                        f"{seg.get('end', seg.get('loc.end',''))}<br>"
                        f"Class: {cls}<br>TCN: {tcn}  LCN: {lcn}"
                    ),
                    hoverinfo="text",
                    showlegend=False,
                ),
                row=row_idx, col=1,
            )

    add_cnv_traces(pri, 1)
    add_cnv_traces(rel, 2)

    # Add chromosome boundary lines
    chr_starts = sorted(
        [(c, v) for c, v in offsets.items()],
        key=lambda x: _chr_to_int(x[0]),
    )
    for chrom, offset in chr_starts:
        for row_idx in [1, 2]:
            fig.add_vline(x=offset, line_width=0.5, line_color="black",
                          row=row_idx, col=1)

    # Chromosome label ticks
    tick_vals  = []
    tick_texts = []
    chr_max_map = {}
    for _, row in all_segs.iterrows():
        c   = str(row["chrom"])
        end = int(row.get("end", row.get("loc.end", 0)))
        chr_max_map[c] = max(chr_max_map.get(c, 0), end)
    for chrom, offset in chr_starts:
        mid = offset + chr_max_map.get(chrom, 0) // 2
        tick_vals.append(mid)
        tick_texts.append(f"chr{chrom}")

    # Legend patches
    for cls, color in CNV_COLORS.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=12, color=color, symbol="square"),
                name=cls,
                showlegend=True,
            ),
            row=1, col=1,
        )

    fig.update_xaxes(
        tickvals=tick_vals,
        ticktext=tick_texts,
        tickangle=45,
        row=2, col=1,
    )
    fig.update_yaxes(visible=False)
    fig.update_layout(
        title=f"{args.patient_id} — Copy Number Alterations",
        height=600,
        template="plotly_white",
        showlegend=True,
    )

    fig.write_html(args.output)
    _try_pdf(fig, args.output)


def _try_pdf(fig, html_path):
    """Export PDF if kaleido is available."""
    try:
        import kaleido  # noqa: F401
        pdf_path = html_path.replace(".html", ".pdf")
        fig.write_image(pdf_path, engine="kaleido")
    except Exception:
        pass


# ---- Entry point ---------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    snv_p = sub.add_parser("snv")
    snv_p.add_argument("--variants",   required=True)
    snv_p.add_argument("--patient-id", required=True)
    snv_p.add_argument("--output",     required=True)

    cnv_p = sub.add_parser("cnv")
    cnv_p.add_argument("--primary-cnv", required=True)
    cnv_p.add_argument("--relapse-cnv", required=True)
    cnv_p.add_argument("--patient-id",  required=True)
    cnv_p.add_argument("--output",      required=True)

    args = p.parse_args()

    if args.command == "snv":
        plot_snv(args)
    elif args.command == "cnv":
        plot_cnv(args)


if __name__ == "__main__":
    main()
