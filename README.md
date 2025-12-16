# RNA-seq analysis – Obeso1 vs Obeso2 (UNIR)

This repository contains a complete RNA-seq workflow using simulated data provided by UNIR.

## Study design
Comparison between two conditions:
- Obeso1: AbrahamSimpson, HomerSimpson
- Obeso2: MargeSimpson, PattyBouvier, SelmaBouvier

## Workflow
1. Quality control of raw reads using FastQC and MultiQC
2. Pseudoalignment and quantification using Salmon
3. Differential expression analysis
4. Biological interpretation and poster

## Repository structure
- `Fastqs/`: raw FASTQ files
- `QC/Raw/`: FastQC and MultiQC reports
- `Salmon/`: Salmon quantification results
- `R/`: R scripts for downstream analysis
- `Figures/`: plots and figures
- `Poster/`: final scientific poster
