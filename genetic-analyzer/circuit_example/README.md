# Genetic circuit RNA-seq analysis — Nextflow DSL2 pipeline

This repository contains a **Nextflow DSL2** workflow to characterize and debug genetic circuits using RNA-seq data. It orchestrates read mapping, feature counting, cohort-level normalization, per-sample transcription profiles, and cryptic site detection using the analysis code in `bin/`.

## Overview

- **Entry point:** `main.nf`
- **Example layout:** `circuit_example/` (contains `data/settings.txt`, `results/`, `tmp/`)
- **Settings file:** tab-delimited `data/settings.txt` (see format below)

### Pipeline stages

1. **Setup** – create `results/` and `tmp/` per sample  
2. **MapReads** – align reads via `bin/map_reads.py` (uses BWA/Samtools)  
3. **CountReads** – per-sample counts, mapped reads, gene lengths via `bin/count_reads.py`  
4. **FragmentDist** – per-sample fragment length distributions via `bin/fragment_distributions.py`  
5. **ReadAnalysis (cohort)** – runs **once** to collate all samples, compute normalization factors, and build matrices via `bin/read_analysis.py` + `binnorm_fpkm.r`  
6. **DistributeNormFactors** – copies shared norm factors to each sample’s results folder  
7. **TranscriptionProfile** – raw and normalized transcription profiles via `bin/transcription_profile.py`  
8. **CrypticSites (cohort)** – runs **once** after all normalized profiles are ready; arranges expected layout and calls `bin/detect_cryptic_tss_tes.py` via the wrapper `bin/cryptic_sites_nf.py` to:
   - detect **cryptic promoters/terminators**,
   - plot selected windows with forward/reverse profiles and feature “rugs”,
   - write a **single** TSV combining known (GFF) sites with detected cryptic sites.

## Requirements

- **Nextflow** ≥ 25.04.x  
- **Java** 8+  
- **Python** 3.9+ with packages used by `bin/*.py`, plus for cryptic detection:
  - `numpy`, `pandas`, `matplotlib`
- **R** 4.x with **edgeR**  
- **BWA** 0.7.x, **Samtools** 1.x, **HTSeq** 0.9.x  

**Optional (macOS quality-of-life):**
- **Graphviz** for DAGs: `brew install graphviz`  
- **Docker Desktop** if you want CPU/memory metrics in `trace.txt` (Nextflow can’t collect them natively on macOS)

## Input data & `settings.txt` format

Place inputs under `data/` and configure a **tab-delimited** file `data/settings.txt`. The first row is the header; the second row is a `None` row that sets a default results directory. Each subsequent row describes a sample.

**Columns (tab-separated):**


**Important:** ensure true **tab** separation (not spaces), and clean sample tokens (e.g., `tube_2`, no trailing spaces).

## Quick start

From the pipeline directory (e.g., `circuit_example/`):

```bash
# optional: clean previous cache/work
nextflow clean -f

# full run with reports
nextflow run main.nf \
  --settings ./data/settings.txt \
  --bin_path ../bin \
  -with-report   report.html \
  -with-timeline timeline.html \
  -with-trace    trace.txt \
  -with-dag      flowchart.svg

# resuming workflow
nextflow run main.nf \
  --settings ./data/settings.txt \
  --bin_path ../bin \
  -resume
