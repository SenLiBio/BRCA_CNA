# CNV Pre-processing Pipeline

---

## Overview

This pipeline performs pre-processing steps for copy number variation (CNV) analysis.  
It converts raw paired-end FASTQ files into cleaned and standardized BAM files suitable for downstream CNV analysis using tools such as **CNVkit** and **Copykit**.

---

## Pipeline Summary

- **Input:** Raw paired-end FASTQ files and the human reference genome (GRCh38)  
- **Steps:**
  1. Read quality control and adapter trimming
  2. Sequence alignment to the reference genome
  3. PCR duplicate marking and removal  
- **Output:** Sorted and deduplicated BAM files ready for input to CNVkit or Copykit

---

## Software Requirements

- `fastp` v0.23.4  
- `bwa` v0.7.17  
- `samtools` v1.19  
- `GATK` v4.3.0.0  

---

## Example Commands

> **Note:** The following example uses `P1.R1.fastq` and `P1.R2.fastq` as input paired-end FASTQ files.

```bash
# Quality control and adapter trimming
fastp \
  -i P1.R1.fastq -I P1.R2.fastq \
  -o P1.R1.clean.fastq -O P1.R2.clean.fastq \
  -j P1.json \
  --detect_adapter_for_pe \
  --cut_front --cut_tail --cut_mean_quality 30

# Alignment to reference genome
bwa mem REFERENCE.fa P1.R1.clean.fastq P1.R2.clean.fastq > P1.sam

# Conversion to BAM and sorting
samtools view -Sb P1.sam | samtools sort -o P1.sorted.bam

# Duplicate marking and indexing
gatk MarkDuplicates \
  --INPUT P1.sorted.bam \
  --OUTPUT P1.dedup.bam \
  --METRICS_FILE P1.dedup.metrics \
  --CREATE_INDEX true
