# CNV Analysis and Visualization for "xxx"

## Overview

This section documents the CNV calling and visualization process using [CNVkit](https://cnvkit.readthedocs.io/en/stable/).  
All analyses are based on BAM files generated in the pre-processing step. Visualizations and downstream analyses were performed using CNVkit’s built-in tools.

---

## Analysis

### Reference and Annotation

- Reference genome: **GRCh38**
- Gene annotation file: `refFlat.txt` (downloaded from UCSC)
- Excluded regions: `excludes.bed` (downloaded from UCSC)

### Example Sample

All commands below use `P1.rmdup.bam` as an example pre-processed BAM file.

### CNVkit Workflow

```bash
# 1. Define accessible regions (excluding known problematic areas)
cnvkit.py access REFERENCE.fa -x excludes.bed -o access-excludes.hg38.bed

# 2. Estimate appropriate bin sizes for WGS data
cnvkit.py autobin P1.rmdup.bam -m wgs -g access-excludes.hg38.bed --annotate refFlat.txt
```

**Note:**  
Based on the estimated resolution from `autobin`, samples were grouped and assigned a fixed bin size as follows:

| Estimated Resolution      | Assigned Bin Size |
|---------------------------|-------------------|
| ≤ 50,000 bp               | 50 kb             |
| 50,000 – 100,000 bp       | 100 kb            |
| 100,000 – 220,000 bp      | 220 kb            |
| 220,000 – 500,000 bp      | 500 kb            |
| 500,000 – 1,000,000 bp    | 1 MB              |
| > 1 MB                    | *Excluded due to excessive noise* |

Example analysis below uses a bin size of **50 kb**:

```bash
# 3. Generate fixed-size bin regions
cnvkit.py autobin P1.rmdup.bam -m wgs \
  -g access-excludes.hg38.bed \
  --annotate refFlat.txt \
  --target-max-size 50000 \
  --target-min-size 50000 \
  --short-names

# 4. Calculate coverage in the target and antitarget regions
cnvkit.py coverage P1.rmdup.bam P1.target.bed -o P1.targetcoverage.cnn
cnvkit.py coverage P1.rmdup.bam P1.antitarget.bed -o P1.antitargetcoverage.cnn

# 5. Build reference from pooled normal samples
cnvkit.py reference *normal.coverage.cnn -f REFERENCE.fa -o Reference.cnn

# 6. Correct raw coverage values for systematic biases
cnvkit.py fix P1.targetcoverage.cnn P1.antitargetcoverage.cnn Reference.cnn -o P1.cnr

# 7. Perform segmentation
cnvkit.py segment P1.cnr --drop-low-coverage -m cbs --smooth-cbs -o P1.cns

# 8. Compute segment-level metrics
cnvkit.py segmetrics P1.cnr -s P1.cns --sem -o P1.segmetrics.cns

# 9. Call absolute copy numbers (with filtering by standard error)
cnvkit.py call --filter sem --drop-low-coverage P1.segmetrics.cns -o P1.call.sem.cns
```

---

