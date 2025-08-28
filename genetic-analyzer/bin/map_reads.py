#!/usr/bin/env python3

import argparse, os, subprocess
import pandas as pd

def parse_settings(settings_path):
    return pd.read_csv(settings_path, sep="\t", index_col=0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-settings", required=True)
    parser.add_argument("-samples", required=True)
    args = parser.parse_args()

    settings = parse_settings(args.settings)
    sample = args.samples
    s = settings.loc[sample]

    fastq1 = s["R1_fastq_file"]
    fastq2 = s.get("R2_fastq_file", "")
    fasta = s["fasta_file"]
    gff   = s["gff_file"]

    ref_base = os.path.splitext(os.path.basename(fasta))[0]
    tmp_dir = f"tmp/{sample}"
    os.makedirs(tmp_dir, exist_ok=True)

    ref_prefix = os.path.join(tmp_dir, sample)

    # Build BWA index
    print(f"Making index: bwa index -p {ref_prefix} {fasta}")
    subprocess.run(["bwa", "index", "-p", ref_prefix, fasta], check=True)

    # Align
    sam_path = os.path.join(tmp_dir, f"{sample}.sam")
    if fastq2 == "" or pd.isna(fastq2):
        bwa_cmd = ["bwa", "mem", ref_prefix, fastq1]
    else:
        bwa_cmd = ["bwa", "mem", ref_prefix, fastq1, fastq2]

    with open(sam_path, "w") as sam:
        print("Running alignment:", " ".join(bwa_cmd))
        subprocess.run(bwa_cmd, stdout=sam, check=True)

    # Convert to BAM
    bam_path = os.path.join(tmp_dir, f"{sample}.bam")
    subprocess.run(["samtools", "view", "-bS", sam_path, "-o", bam_path], check=True)
    subprocess.run(["samtools", "sort", bam_path, "-o", bam_path], check=True)
    subprocess.run(["samtools", "index", bam_path], check=True)

if __name__ == "__main__":
    main()

