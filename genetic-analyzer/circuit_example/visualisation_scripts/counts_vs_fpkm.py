#!/usr/bin/env python3
"""
Scatter plot combining counts and FPKM for each gene (per sample).

X-axis: raw counts (from *.counts.txt)
Y-axis: log2(FPKM+1) (from fpkm.normed.matrix.txt)

Output:
  results/<sample>/<sample>_counts_vs_fpkm.png
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MATRIX = "results/fpkm.normed.matrix.txt"

def main():
    fpkm = pd.read_csv(MATRIX, sep="\t", index_col=0)
    fpkm_log = np.log2(fpkm.clip(lower=0) + 1)

    for counts_file in glob.glob("results/*/*.counts.txt"):
        sample = os.path.basename(counts_file).replace(".counts.txt", "")
        outdir = os.path.join("results", sample)
        os.makedirs(outdir, exist_ok=True)

        # Load counts
        counts = pd.read_csv(counts_file, sep="\t", header=None,
                             names=["gene", "count"], dtype={"gene": str})
        counts = counts[~counts["gene"].str.startswith("__")]
        counts.set_index("gene", inplace=True)

        if sample not in fpkm.columns:
            print(f"[WARN] {sample} not in FPKM matrix, skipping.")
            continue

        # Join counts + FPKM
        merged = counts.join(fpkm_log[[sample]], how="inner")
        merged.columns = ["counts", "log2(FPKM+1)"]

        # Scatter plot
        plt.figure(figsize=(7, 6))
        plt.scatter(merged["counts"], merged["log2(FPKM+1)"],
                    alpha=0.4, s=10)
        plt.xscale("log")
        plt.xlabel("Raw counts (log scale)")
        plt.ylabel("log2(FPKM+1)")
        plt.title(f"Counts vs FPKM — {sample}")
        plt.tight_layout()

        out_png = os.path.join(outdir, f"{sample}_counts_vs_fpkm.png")
        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()

