#!/usr/bin/env python3
"""
cryptic_sites_nf.py — Nextflow-friendly wrapper for detect_cryptic_tss_tes.py

What it expects (in the task work dir):
  *.fwd.norm.profiles.txt
  *.rev.norm.profiles.txt

What it does:
  - Infers sample IDs from filenames (prefix before the first '.').
  - Creates results/<sample>/ and copies those norm-profile files there.
  - Calls detect_cryptic_tss_tes.py with --base-dir .
  - Leaves 'all_promoters_terminators_with_cryptics.txt' in CWD.
  - Plots (PNGs) are emitted under results/<sample>/ by your script.
"""

import argparse, sys, shutil, subprocess
from pathlib import Path
from collections import defaultdict


def stage_profiles_into_results():
    cwd = Path.cwd()
    profiles = list(cwd.glob("*.norm.profiles.txt"))
    if not profiles:
        print("ERROR: No *.norm.profiles.txt files found in the working directory.", file=sys.stderr)
        sys.exit(2)

    by_sample = defaultdict(list)
    for p in profiles:
        sample = p.name.split(".", 1)[0]
        by_sample[sample].append(p)

    staged = []
    for sample, files in by_sample.items():
        dest_dir = cwd / "results" / sample
        dest_dir.mkdir(parents=True, exist_ok=True)
        for f in files:
            dest = dest_dir / f.name
            shutil.copy2(str(f), str(dest))
            staged.append(dest)
    return staged


def find_detector():
    here = Path(__file__).resolve().parent
    detector = here / "detect_cryptic_tss_tes.py"
    if not detector.exists():
        print(f"ERROR: Cannot find detector script at {detector}", file=sys.stderr)
        sys.exit(3)
    return str(detector)


def main():
    ap = argparse.ArgumentParser(description="Nextflow wrapper for cryptic TSS/TES detection")
    ap.add_argument("--gff", required=True, help="Path to genome annotations (GFF/GFF3/GTF)")
    args = ap.parse_args()

    staged = stage_profiles_into_results()
    print(f"Staged {len(staged)} profile files under results/<sample>/")

    detector = find_detector()
    cmd = [sys.executable, detector, "--base-dir", ".", "--gff", args.gff]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    combined = Path("all_promoters_terminators_with_cryptics.txt")
    if combined.exists():
        print("Wrote:", combined.resolve())
    else:
        print("WARNING: Combined TSV not found; detector may have written a different name.", file=sys.stderr)


if __name__ == "__main__":
    main()

