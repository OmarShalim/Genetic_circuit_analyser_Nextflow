#!/usr/bin/env python3
"""
Plot top-N gene counts for all *.counts.txt under results/.

- Searches:
    results/*.counts.txt
    results/*/*.counts.txt
- Filters out technical categories starting with "__" (e.g., __ambiguous).
- Saves:
    results/<sample>/<sample>_top{N}_counts.png

Usage:
  python count_reads_top20_plot.py            # default top_n=20
  python count_reads_top20_plot.py 30         # top 30
"""

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt


def plot_top_counts(infile: str, top_n: int = 20) -> str:
    # sample name from filename without extension
    base = os.path.basename(infile)
    sample = base.replace(".counts.txt", "")
    outdir = os.path.join("results", sample)
    os.makedirs(outdir, exist_ok=True)

    # Read counts (gene \t count). Be tolerant of odd rows.
    df = pd.read_csv(infile, sep="\t", header=None, names=["gene", "count"], dtype={"gene": str}, engine="python")

    # Drop technical/non-gene categories (HTSeq/featureCounts) that start with "__"
    df = df[~df["gene"].astype(str).str.startswith("__")]

    # Keep only numeric counts
    df["count"] = pd.to_numeric(df["count"], errors="coerce")
    df = df.dropna(subset=["count"])

    if df.empty:
        print(f"[WARN] No biological genes found in {infile} after filtering.")
        return ""

    # Top N
    df_top = df.sort_values("count", ascending=False).head(top_n)

    # Plot
    plt.figure(figsize=(12, 7))
    plt.barh(df_top["gene"][::-1], df_top["count"][::-1])
    plt.xlabel("Counts")
    plt.ylabel("Gene")
    plt.title(f"Top {len(df_top)} counts - {sample}")
    plt.tight_layout()

    out_png = os.path.join(outdir, f"{sample}_top{len(df_top)}_counts.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    return out_png


def main():
    # Optional CLI arg: top_n
    top_n = 20
    if len(sys.argv) > 1:
        try:
            top_n = int(sys.argv[1])
        except ValueError:
            pass

    # Find files at results/*.counts.txt and results/*/*.counts.txt
    files = sorted(set(glob.glob("results/*.counts.txt") + glob.glob("results/*/*.counts.txt")))
    if not files:
        print("No *.counts.txt files found under results/ or results/*/")
        return

    for infile in files:
        out_png = plot_top_counts(infile, top_n=top_n)
        if out_png:
            print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()

