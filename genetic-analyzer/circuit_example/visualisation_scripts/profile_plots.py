#!/usr/bin/env python3
"""
Plot forward and reverse strand profiles for each sample.

Looks for (one level below results/):
  results/<sample>/<sample>.fwd.norm.profiles.txt
  results/<sample>/<sample>.rev.norm.profiles.txt

Assumed TSV columns (no header):
  contig   start   end   position   value

Output:
  results/<sample>/<sample>_strand_profiles.png

Default: smoothing window = 31
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt


def read_profile(path: str) -> pd.DataFrame:
    """Read a strand profile TSV without header and return positions+values (sorted)."""
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["contig", "start", "end", "pos", "value"]
    )
    df = df[["pos", "value"]].apply(pd.to_numeric, errors="coerce").dropna()
    df = df.groupby("pos", as_index=False)["value"].mean().sort_values("pos")
    return df


def smooth_series(y: pd.Series, window: int = 31) -> pd.Series:
    """Centered rolling mean (default window=31)."""
    return y.rolling(window=window, center=True, min_periods=1).mean()


def plot_sample(sample_dir: str, sample: str, window: int = 31) -> str:
    fwd_path = os.path.join(sample_dir, f"{sample}.fwd.norm.profiles.txt")
    rev_path = os.path.join(sample_dir, f"{sample}.rev.norm.profiles.txt")

    has_fwd = os.path.exists(fwd_path)
    has_rev = os.path.exists(rev_path)
    if not (has_fwd or has_rev):
        return ""

    os.makedirs(sample_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, figsize=(12, 6),
        gridspec_kw=dict(hspace=0.1)
    )

    if has_fwd:
        df_f = read_profile(fwd_path)
        ax1.plot(df_f["pos"], smooth_series(df_f["value"], window), color="tab:blue")
        ax1.set_ylabel("Forward")
        ax1.set_title(f"Strand profiles — {sample}")

    if has_rev:
        df_r = read_profile(rev_path)
        ax2.plot(df_r["pos"], smooth_series(df_r["value"], window), color="tab:orange")
        ax2.set_ylabel("Reverse")
        ax2.set_xlabel("Position")

    out_png = os.path.join(sample_dir, f"{sample}_strand_profiles.png")
    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return out_png


def main():
    sample_dirs = [d for d in glob.glob("results/*") if os.path.isdir(d)]
    if not sample_dirs:
        print("No sample directories found under results/*")
        return

    for d in sorted(sample_dirs):
        sample = os.path.basename(d)
        out_png = plot_sample(d, sample, window=31)
        if out_png:
            print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
