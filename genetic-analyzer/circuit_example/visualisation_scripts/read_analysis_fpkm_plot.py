#!/usr/bin/env python3
"""
FPKM visuals:
  - Per-sample histogram of log2(FPKM+1) (optional KDE overlay)
  - Combined violin plot across samples

Inputs:
  results/fpkm.normed.matrix.txt  (rows = genes, cols = samples)

Outputs:
  results/<sample>/<sample>_fpkm_hist.png
  results/fpkm_violin.png
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Violin needs seaborn
import seaborn as sns

# Optional KDE for histograms
try:
    from scipy.stats import gaussian_kde
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def plot_sample_hist(sample: str, values: pd.Series, outdir: str, bins: int = 60, kde: bool = False) -> str:
    """
    Plot histogram of log2(FPKM+1) for one sample and save PNG.
    Y-axis = number of genes (counts). Optional KDE overlay (scaled to counts).
    """
    log_vals = np.log2(values.clip(lower=0) + 1)
    log_vals = log_vals.replace([np.inf, -np.inf], np.nan).dropna()

    plt.figure(figsize=(8, 6))
    counts, edges, _ = plt.hist(log_vals, bins=bins, alpha=0.85)
    plt.xlabel("log2(FPKM + 1)")
    plt.ylabel("Number of genes")
    plt.title(f"Expression distribution — {sample}")

    if kde and _HAVE_SCIPY and len(log_vals) > 1:
        kde_model = gaussian_kde(log_vals)
        x_grid = np.linspace(edges[0], edges[-1], 512)
        bin_width = (edges[-1] - edges[0]) / bins
        density_counts = kde_model(x_grid) * len(log_vals) * bin_width
        plt.plot(x_grid, density_counts, linewidth=2)

    plt.tight_layout()
    out_png = os.path.join(outdir, f"{sample}_fpkm_hist.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    return out_png


def plot_violin(df: pd.DataFrame, out_png: str) -> str:
    """Make a combined violin plot across samples for log2(FPKM+1)."""
    df_log = np.log2(df.clip(lower=0) + 1)
    df_long = df_log.melt(var_name="Sample", value_name="log2(FPKM+1)")
    df_long = df_long.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2(FPKM+1)"])

    plt.figure(figsize=(max(12, 1.2 * len(df.columns)), 6))
    sns.violinplot(
        data=df_long,
        x="Sample",
        y="log2(FPKM+1)",
        inner="quartile",
        scale="width",
        cut=0
    )
    plt.title("FPKM distributions across samples")
    plt.ylabel("log2(FPKM + 1)")
    plt.xlabel("Samples")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    ensure_dir(os.path.dirname(out_png) or ".")
    plt.savefig(out_png, dpi=300)
    plt.close()
    return out_png


def main():
    parser = argparse.ArgumentParser(description="FPKM histogram(s) and violin plot")
    parser.add_argument("--matrix", default="results/fpkm.normed.matrix.txt",
                        help="Path to fpkm.normed.matrix.txt (default: results/fpkm.normed.matrix.txt)")
    parser.add_argument("--bins", type=int, default=60, help="Histogram bins (default: 60)")
    parser.add_argument("--kde", action="store_true", help="Overlay KDE on histograms (SciPy required)")
    parser.add_argument("--skip-hist", action="store_true", help="Skip per-sample histograms")
    parser.add_argument("--skip-violin", action="store_true", help="Skip combined violin plot")
    args = parser.parse_args()

    # Load matrix (genes x samples)
    df = pd.read_csv(args.matrix, sep="\t", index_col=0)

    # Histograms per sample
    if not args.skip_hist:
        for sample in df.columns:
            outdir = os.path.join("results", sample)
            ensure_dir(outdir)
            out_png = plot_sample_hist(sample, df[sample], outdir, bins=args.bins, kde=args.kde)
            print(f"Saved: {out_png}")
        if args.kde and not _HAVE_SCIPY:
            print("Note: SciPy not found; KDE overlay was skipped. Install with: pip install scipy")

    # Combined violin
    if not args.skip_violin:
        out_png = os.path.join("results", "fpkm_violin.png")
        out_png = plot_violin(df, out_png)
        print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
