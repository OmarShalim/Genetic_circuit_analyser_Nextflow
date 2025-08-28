#!/usr/bin/env python3
"""
Detect cryptic TSS/TES, plot forward/reverse profiles, and write ONE TSV
combining known (GFF) promoters/terminators with detected cryptic sites.

Directory layout supported (example):
  .
  ├─ results/
  │   ├─ tube_1/  *.fwd.norm.profiles.txt, *.rev.norm.profiles.txt
  │   ├─ tube_2/  ...
  │   └─ tube_N/  ...
  └─ data/gff/   your_annotations.gff

Run from your project root (or pass --base-dir).
Plots are written into each tube folder; the combined TSV is written to --base-dir.

Outputs:
  - results/tube_X/<sample>_<seqid>_<start>-<end>_clean_rugs_colored.png  (3 windows per tube)
  - <base-dir>/all_promoters_terminators_with_cryptics.txt  (TSV; known + cryptic)

Usage:
  python detect_cryptic_tss_tes.py
  python detect_cryptic_tss_tes.py --base-dir /path/to/project --gff /path/to/project/data/gff/anno.gff
  python detect_cryptic_tss_tes.py --windows 4 --margin 2000
"""

import os, glob, re, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------- Appearance & detection defaults ----------------
COLOR_FWD  = "tab:orange"
COLOR_REV  = "tab:cyan"
COLOR_PROM = "tab:blue"
COLOR_TERM = "tab:red"
COLOR_CPR  = "tab:green"   # cryptic promoter
COLOR_CTR  = "tab:purple"  # cryptic terminator

SMOOTH_DISPLAY_WIN = 21           # bp for plotting only
DETECT_SMOOTH_WIN  = 31           # bp for cryptic detection
K_START, K_MIN, K_STEP = 3.0, 0.6, 0.3
EXCLUDE_KNOWN_BP   = 120          # don't call cryptic too close to known

FWD_SUF = ".fwd.norm.profiles.txt"
REV_SUF = ".rev.norm.profiles.txt"
TUBE_DIR_REGEX = re.compile(r"^tube_\d+$")

# ---------------- File I/O helpers ----------------
def read_profile(path: str) -> pd.DataFrame:
    """Read a norm profile: seqid, region_start, region_end, position, value."""
    cols = ["seqid","region_start","region_end","position","value"]
    df = pd.read_csv(path, sep=r"\s+|,", engine="python", header=None, names=cols)
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["value"]    = pd.to_numeric(df["value"],    errors="coerce")
    df = df.dropna(subset=["position","value"]).reset_index(drop=True)
    df["position"] = df["position"].astype(int)
    return df

def parse_gff(path: str) -> pd.DataFrame:
    """Parse GFF and extract promoter/terminator features with centers."""
    cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]
    rows = []
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=cols+["kind","center"])
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                parts += [""]*(9-len(parts))
            row = dict(zip(cols, parts[:9]))
            try:
                row["start"] = int(float(row["start"]))
                row["end"]   = int(float(row["end"]))
            except Exception:
                continue
            rows.append(row)
    gff = pd.DataFrame(rows)

    def classify(t, a):
        t = (t or "").lower()
        a = (a or "").lower()
        if "promoter" in t or "promoter" in a:
            return "promoter"
        if "terminator" in t or "terminator" in a or "rho_independent_terminator" in a or "rho-independent" in a:
            return "terminator"
        return None

    gff["kind"] = [classify(t, a) for t, a in zip(gff["type"], gff["attributes"])]
    gff = gff[gff["kind"].notna()].copy()
    if not gff.empty:
        gff["center"] = (gff["start"] + gff["end"]) // 2
    return gff.sort_values(["seqid","start","end"]).reset_index(drop=True)

# ---------------- Recursive pairing over results/tube_* ----------------
def _sample_key(path: str) -> str:
    """Prefer the parent dir if it looks like 'tube_X'; else use filename prefix."""
    parent = os.path.basename(os.path.dirname(path))
    base   = os.path.basename(path)
    if TUBE_DIR_REGEX.match(parent):
        return parent
    if base.endswith(FWD_SUF):
        return base[:-len(FWD_SUF)]
    if base.endswith(REV_SUF):
        return base[:-len(REV_SUF)]
    return parent  # fallback

def find_pairs(base_dir: str):
    """
    Recursively find fwd/rev files and pair them by tube folder or filename prefix.
    Returns list of tuples: (key, fwd_path, rev_path)
    """
    fwd_globs = glob.glob(os.path.join(base_dir, "**", f"*{FWD_SUF}"), recursive=True)
    rev_globs = glob.glob(os.path.join(base_dir, "**", f"*{REV_SUF}"), recursive=True)

    fwd_map, rev_map = {}, {}
    for p in fwd_globs:
        key = _sample_key(p)
        # prefer a file inside a tube_* folder if duplicates
        if key not in fwd_map or TUBE_DIR_REGEX.match(os.path.basename(os.path.dirname(p))):
            fwd_map[key] = p
    for p in rev_globs:
        key = _sample_key(p)
        if key not in rev_map or TUBE_DIR_REGEX.match(os.path.basename(os.path.dirname(p))):
            rev_map[key] = p

    # Additionally ensure tube_* directories are considered as keys
    for p in fwd_globs + rev_globs:
        parent = os.path.basename(os.path.dirname(p))
        if TUBE_DIR_REGEX.match(parent):
            if parent not in fwd_map:
                cand = next((x for x in fwd_globs if os.path.basename(os.path.dirname(x)) == parent), None)
                if cand: fwd_map[parent] = cand
            if parent not in rev_map:
                cand = next((x for x in rev_globs if os.path.basename(os.path.dirname(x)) == parent), None)
                if cand: rev_map[parent] = cand

    pairs = []
    for k in sorted(set(fwd_map) & set(rev_map)):
        if fwd_map[k] and rev_map[k]:
            pairs.append((k, fwd_map[k], rev_map[k]))
    return pairs

# ---------------- Detection ----------------
def moving_average(y, win=31):
    win = max(3, int(win) | 1)  # odd, >=3
    s = pd.Series(y)
    return s.rolling(win, center=True, min_periods=1).mean().values

def detect_cryptic_mad(pos, fwd, rev, known_proms, known_terms,
                       exclude_bp=EXCLUDE_KNOWN_BP, smooth_win=DETECT_SMOOTH_WIN,
                       k_start=K_START, k_min=K_MIN, k_step=K_STEP) -> pd.DataFrame:
    """
    Find cryptic sites as large +/- gradients in smoothed (fwd+rev).
    Returns DataFrame with columns: type, position, score
    """
    pos = np.asarray(pos)
    t   = np.asarray(fwd) + np.asarray(rev)
    t_s = moving_average(t, win=smooth_win)
    g   = np.gradient(t_s, pos)

    g_med = np.median(g)
    mad   = np.median(np.abs(g - g_med))
    sigma = 1.4826 * mad if mad > 0 else (np.std(g) if np.std(g) > 0 else 1.0)

    known_centers = np.array([(a+b)//2 for a, b in (known_proms + known_terms)]) \
                    if (known_proms or known_terms) else np.array([], dtype=int)

    def too_close(p):
        return known_centers.size > 0 and np.any(np.abs(known_centers - p) <= exclude_bp)

    k = k_start
    picks_prom, picks_term = [], []
    while k >= k_min and (len(picks_prom) + len(picks_term) < 2):
        pos_thr = g_med + k * sigma
        neg_thr = g_med - k * sigma
        idx_prom = np.where(g >= pos_thr)[0]
        idx_term = np.where(g <= neg_thr)[0]

        def thin(idxs, min_sep=120):
            if len(idxs) == 0:
                return []
            idxs = sorted(idxs, key=lambda i: abs(g[i]), reverse=True)
            keep = []
            for i in idxs:
                p = pos[i]
                if (not too_close(p)) and all(abs(p - pos[j]) >= min_sep for j in keep):
                    keep.append(i)
                if len(keep) >= 200:
                    break
            return keep

        picks_prom = thin(idx_prom)
        picks_term = thin(idx_term)
        if len(picks_prom) + len(picks_term) >= 2:
            break
        k -= k_step

    rows = []
    for i in picks_prom:
        rows.append({"type":"cryptic_promoter",  "position": int(pos[i]), "score": float(g[i])})
    for i in picks_term:
        rows.append({"type":"cryptic_terminator","position": int(pos[i]), "score": float(g[i])})
    rows.sort(key=lambda d: abs(d["score"]), reverse=True)
    return pd.DataFrame(rows)

# ---------------- Window selection & plotting ----------------
def choose_windows(df_fwd, df_rev, cryptic_df, n=3, margin=1500):
    """Choose up to n windows around cryptic clusters; fallback to strongest peaks."""
    windows = []
    if cryptic_df is not None and not cryptic_df.empty:
        ps = sorted(cryptic_df["position"].astype(int).tolist())
        clusters = [[ps[0]]]
        for p in ps[1:]:
            if p - clusters[-1][-1] <= 1500:
                clusters[-1].append(p)
            else:
                clusters.append([p])
        clusters = sorted(clusters, key=lambda xs: (len(xs), xs[-1] - xs[0]), reverse=True)
        pmin, pmax = int(df_fwd["position"].min()), int(df_fwd["position"].max())
        for xs in clusters[:n]:
            windows.append({"start": max(pmin, xs[0]-margin), "end": min(pmax, xs[-1]+margin)})

    if len(windows) < n:
        merged = df_fwd[["position","value"]].merge(
            df_rev[["position","value"]], on="position", suffixes=("_fwd","_rev"))
        merged["sum"] = merged["value_fwd"] + merged["value_rev"]
        merged = merged.sort_values("sum", ascending=False)
        seen = []
        for _, row in merged.iterrows():
            p = int(row["position"])
            if any((p >= w["start"] and p <= w["end"]) for w in windows):
                continue
            if all(abs(p - q) >= 3000 for q in seen):
                windows.append({"start": max(int(df_fwd["position"].min()), p-margin),
                                "end":   min(int(df_fwd["position"].max()), p+margin)})
                seen.append(p)
            if len(windows) == n:
                break
    return windows[:n]

def plot_window_clean_rugs_colored(prefix, seqid, df_fwd, df_rev, gff, cryptic_df,
                                   window, outdir=".", smooth_win=SMOOTH_DISPLAY_WIN):
    """Draw clean profiles with feature 'rug' ticks in a small band below the plot."""
    wstart, wend = int(window["start"]), int(window["end"])
    fwd = df_fwd[(df_fwd["position"] >= wstart) & (df_fwd["position"] <= wend)].copy()
    rev = df_rev[(df_rev["position"] >= wstart) & (df_rev["position"] <= wend)].copy()
    if fwd.empty or rev.empty:
        return None

    fwd_s = moving_average(fwd["value"].values, win=smooth_win)
    rev_s = moving_average(rev["value"].values, win=smooth_win)

    fig, ax = plt.subplots(figsize=(11.5, 5))
    ax.plot(fwd["position"], fwd_s, linewidth=1.2, alpha=0.95, color=COLOR_FWD, label="Forward profile")
    ax.plot(rev["position"], rev_s, linewidth=1.2, alpha=0.95, color=COLOR_REV, label="Reverse profile")

    # Reserve a bottom band for rug ticks
    ymin, ymax = ax.get_ylim()
    yr = ymax - ymin
    band_height = 0.05 * yr
    gap         = 0.01 * yr
    n_tracks    = 4
    extra_space = n_tracks*(band_height+gap) + gap
    ax.set_ylim(ymin - extra_space, ymax)
    base = ymin - extra_space + gap

    def rug(xs, track_idx, color="k", lw=1.6):
        y0 = base + track_idx*(band_height+gap)
        for x in xs:
            ax.vlines(x, y0, y0+band_height, color=color, linewidth=lw)

    # Known features as rug ticks (centers)
    gf_sub = gff[(gff["seqid"] == seqid) & (gff["start"] <= wend) & (gff["end"] >= wstart)]
    prom_centers = [int((int(r["start"])+int(r["end"]))//2) for _, r in gf_sub[gf_sub["kind"]=="promoter"].iterrows()]
    term_centers = [int((int(r["start"])+int(r["end"]))//2) for _, r in gf_sub[gf_sub["kind"]=="terminator"].iterrows()]
    rug(prom_centers, track_idx=0, color=COLOR_PROM)
    rug(term_centers, track_idx=1, color=COLOR_TERM)

    # Cryptic ticks
    if cryptic_df is not None and not cryptic_df.empty:
        cry = cryptic_df[(cryptic_df["position"]>=wstart) & (cryptic_df["position"]<=wend)]
        rug(cry.loc[cry["type"]=="cryptic_promoter","position"].astype(int).tolist(),   track_idx=2, color=COLOR_CPR)
        rug(cry.loc[cry["type"]=="cryptic_terminator","position"].astype(int).tolist(), track_idx=3, color=COLOR_CTR)

    ax.set_title(f"{prefix}: {wstart}-{wend}")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Normalized signal")

    legend_elements = [
        Line2D([0],[0], color=COLOR_FWD, label="Forward profile"),
        Line2D([0],[0], color=COLOR_REV, label="Reverse profile"),
        Line2D([0],[0], color=COLOR_PROM, label="promoter"),
        Line2D([0],[0], color=COLOR_TERM, label="terminator"),
        Line2D([0],[0], color=COLOR_CPR, label="cryptic promoter"),
        Line2D([0],[0], color=COLOR_CTR, label="cryptic terminator"),
    ]
    fig.subplots_adjust(right=0.82)
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.margins(x=0.02)
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    fname = f"{prefix}_{seqid}_{wstart}-{wend}_clean_rugs_colored.png".replace("/", "_")
    out_path = os.path.join(outdir, fname)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path

# ---------------- Orchestration ----------------
def run_all(base_dir=".", gff_path=None, n_windows=3, margin=1500):
    pairs = find_pairs(base_dir)
    if not pairs:
        print("No fwd/rev pairs found under", base_dir)
        return [], None

    gff = parse_gff(gff_path or os.path.join(base_dir, "0x58v50.gff"))
    image_paths = []
    combined_rows = []

    for key, fwd_path, rev_path in pairs:
        df_fwd = read_profile(fwd_path)
        df_rev = read_profile(rev_path)
        seqid  = df_fwd["seqid"].iloc[0]

        # Known features -> combined_rows
        gf_seq = gff[gff["seqid"]==seqid].copy()
        if not gf_seq.empty:
            gf_seq = gf_seq[["seqid","kind","start","end","strand"]].copy()
            gf_seq["center"] = ((gf_seq["start"] + gf_seq["end"])//2).astype(int)
            for _, r in gf_seq.iterrows():
                combined_rows.append({
                    "category":"known","sample":key,"seqid":r["seqid"],"type":r["kind"],
                    "start":int(r["start"]),"end":int(r["end"]),
                    "center_or_position":int(r["center"]),"strand":str(r.get("strand",".")),"score":""
                })

        # Cryptic detection
        merged = df_fwd[["position","value"]].merge(
            df_rev[["position","value"]], on="position", suffixes=("_fwd","_rev")
        ).sort_values("position")
        pos = merged["position"].values
        fwd_vals = merged["value_fwd"].values
        rev_vals = merged["value_rev"].values
        known_proms = [(int(r["start"]), int(r["end"])) for _, r in gff[(gff["seqid"]==seqid)&(gff["kind"]=="promoter")].iterrows()]
        known_terms = [(int(r["start"]), int(r["end"])) for _, r in gff[(gff["seqid"]==seqid)&(gff["kind"]=="terminator")].iterrows()]
        cryptic_df = detect_cryptic_mad(pos, fwd_vals, rev_vals, known_proms, known_terms)

        for _, c in cryptic_df.iterrows():
            combined_rows.append({
                "category":"cryptic","sample":key,"seqid":seqid,"type":c["type"],
                "start":"", "end":"", "center_or_position":int(c["position"]),
                "strand":".","score":float(c["score"])
            })

        # Windows & plots — save images into the tube's own directory
        outdir = os.path.dirname(fwd_path)
        windows = choose_windows(df_fwd, df_rev, cryptic_df, n=n_windows, margin=margin)
        for w in windows:
            out = plot_window_clean_rugs_colored(key, seqid, df_fwd, df_rev, gff, cryptic_df, w, outdir=outdir)
            if out: image_paths.append(out)

    # Write single TSV at base_dir
    combined_df = pd.DataFrame(combined_rows, columns=[
        "category","sample","seqid","type","start","end","center_or_position","strand","score"
    ]).sort_values(["sample","seqid","category","type","center_or_position"])
    out_txt = os.path.join(base_dir, "all_promoters_terminators_with_cryptics.txt")
    combined_df.to_csv(out_txt, sep="\t", index=False)
    return image_paths, out_txt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", default=".", help="Project root (scanned recursively for results/tube_*/)")
    ap.add_argument("--gff", default=None, help="Path to GFF (default: auto/ <base-dir>/0x58v50.gff)")
    ap.add_argument("--windows", type=int, default=3, help="How many windows to plot per tube")
    ap.add_argument("--margin",  type=int, default=1500, help="Half-width margin around clusters/peaks")
    args = ap.parse_args()

    imgs, txt_path = run_all(base_dir=args.base_dir, gff_path=args.gff, n_windows=args.windows, margin=args.margin)
    print("Saved TXT:", txt_path)
    if imgs:
        print("Saved plots:")
        for p in imgs:
            print(" ", p)

if __name__ == "__main__":
    main()

