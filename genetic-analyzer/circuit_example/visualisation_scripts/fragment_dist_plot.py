#!/usr/bin/env python3
"""
Generate two fragment distribution plots for each sample:

1. 0–1000 bp
2. 100–400 bp

Looks for: results/*/*.fragment.distribution.txt
Saves to:  results/<sample>/<sample>_fragment_distribution_0-1000bp.png
           results/<sample>/<sample>_fragment_distribution_100-400bp.png
"""

import argparse
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt


def get_sample_name(file_path: str) -> str:
    """Use the parent directory (results/<sample>/) as sample name."""
    return os.path.basename(os.path.dirname(file_path))


def read_distribution(file_path: str) -> pd.DataFrame:
    """Read two-column TSV: fragment_length \t count."""
    df = pd.read_csv(file_path, sep="\t", header=None, names=["fragment_length", "count"])
    df = df.dropna()
    df["fragment_length"] = pd.to_numeric(df["fragment_length"], errors="coerce")
    df["count"] = pd.to_numeric(df["count"], errors="coerce")
    df = df.dropna().sort_values("fragment_length")
    return df


def plot_range(df: pd.DataFrame, sample: str, outdir: str, min_bp: int, max_bp: int) -> str:
    """Filter to [min_bp, max_bp] and save a step plot."""
    df_range = df[(df["fragment_length"] >= min_bp) & (df["fragment_length"] <= max_bp)]

    plt.figure(figsize=(10, 6))
    plt.step(df_range["fragment_length"], df_range["count"], where="mid")
    plt.xlabel("Fragment length (bp)")
    plt.ylabel("Read count")
    plt.title(f"Fragment length distribution — {sample} ({min_bp}-{max_bp} bp)")
    plt.tight_layout()

    out_png = os.path.join(outdir, f"{sample}_fragment_distribution_{min_bp}-{max_bp}bp.png")
    plt.savefig(out_png, dpi=300)
    plt.close()
    return out_png


def process_file(file_path: str):
    sample = get_sample_name(file_path)
    outdir = os.path.join("results", sample)
    os.makedirs(outdir, exist_ok=True)

    df = read_distribution(file_path)

    p1 = plot_range(df, sample, outdir, 0, 1000)
    p2 = plot_range(df, sample, outdir, 100, 400)

    print(f"Saved: {p1}")
    print(f"Saved: {p2}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", help="Optional specific .fragment.distribution.txt files to process")
    ns = parser.parse_args()

    files = ns.paths if ns.paths else glob.glob(os.path.join("results", "*", "*.fragment.distribution.txt"))
    if not files:
        print("No .fragment.distribution.txt files found under results/*/.")
        return

    for fp in sorted(files):
        process_file(fp)


if __name__ == "__main__":
    main()
