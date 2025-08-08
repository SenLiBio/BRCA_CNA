# CNV Analysis Using Single-Cell DNA Sequencing Data

## Overview

This repository demonstrates the workflow for calling and visualizing copy number variations (CNVs) from single-cell DNA sequencing data using [CopyKit](https://github.com/navinlabcode/copykit).  
All analyses are based on BAM files generated from pre-processed single-cell data. CNV calling, quality control, clustering, and phylogenetic inference were performed using CopyKit’s built-in functions, along with additional custom R scripts.

> 💡 **Note:** Some custom functions (e.g., `filter_cells_by_sds`) used in this pipeline will be packaged into a standalone R package in the future.

---
## Notes
- This pipeline assumes that the input BAM files are already deduplicated and indexed.
- You can adjust resolution, QC thresholds (SDS ≤ 1, correlation ≥ 0.8), or clustering resolution based on dataset characteristics.
- More details on custom functions and SDS-based QC are provided in the `scripts/` directory.


## Analysis Pipeline

```r
library(copykit)

# 1. Pre-processing
tumor <- runVarbin("~/path/to/scDNA/bam/files/", 
                   remove_Y = TRUE,
                   genome = "hg38",
                   resolution = "220kb",
                   is_paired_end = TRUE,
                   method="multipcf")

# 2. Quality control
tumor <- runMetrics(tumor)  # Add basic QC metrics
tumor <- findAneuploidCells(tumor)  # Identify euploid vs. aneuploid
tumor <- tumor[, colData(tumor)$is_aneuploid == TRUE]

# 3. Outlier filtering (custom + built-in)
tumor <- findOutliers(tumor, k = 5, resolution = 0.8)
tumor <- filter_cells_by_sds(tumor)  # Custom function for SDS-based QC
tumor <- tumor[, colData(tumor)$outlier == FALSE]

# 4. KNN smoothing
tumor <- knnSmooth(tumor)

# 5. Integer copy number inference
tumor <- calcInteger(tumor, method = "scquantum", assay = "smoothed_bincounts")

# 6. Clustering
tumor <- runUmap(tumor)
tumor <- findSuggestedK(tumor)
tumor <- findClusters(tumor)
tumor <- tumor[, colData(tumor)$subclones != "c0"]  # Remove unclustered

# 7. Phylogenetic tree construction
tumor <- runPhylo(tumor, metric = "manhattan")

# 8. Consensus profile generation
tumor <- calcConsensus(tumor)
```

## Visulization
```r
# Add metadata
colData(tumor)$Visit <- stringr::str_extract(colData(tumor)$sample, "(V0|V2)")
colData(tumor)$Info <- stringr::str_extract(colData(tumor)$sample, "\\d{4}_(V0|V2)")
colData(tumor)$Patient <- stringr::str_extract(colData(tumor)$sample, "\\d{4}")

# UMAP plots
plotUmap(tumor, label = "subclones")
plotUmap(tumor, label = "Visit")

# Heatmap: outlier comparison
plotHeatmap(tumor, label = "outlier", row_split = "outlier")

# Heatmap: per-individual consensus profile
P1 <- tumor[, colData(tumor)$sample == "P1"]
plotHeatmap(P1, consensus = TRUE, assay = "integer", label = "subclones", group = "stage")

# Heatmap: full cohort
plotHeatmap(tumor, label = c("Patient", "Visit"), group = "Patient")
```
