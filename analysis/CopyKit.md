# CNV Analysis using scDNA sequencing data

## Overview

This section documents the CNV calling and visualization process using [CopyKit](https://github.com/navinlabcode/copykit?tab=readme-ov-file).  
All analyses are based on BAM files generated in the pre-processing step. Visualizations and downstream analyses were performed using CopyKit’s built-in tools or homemade R scripts.

---

## Analysis

 ***(Later I'll package all the R functions)***
 
```R
library("copykit")

# Pre-processing
tumor <- runVarbin("~/path/to/scDNA/bam/files/", 
                   remove_Y = TRUE,
                   genome = "hg38",
                   resolution = "220kb"
                   is_paired_end = TRUE)
# Quality control
# add basic quality control information to colData
tumor <- runMetrics(tumor)

# detect euploid cells
tumor <- findAneuploidCells(tumor)
tumor <- tumor[,colData(tumor)$is_aneuploid == TRUE]

# annotates low-quality cells according to a defined resolution threshold and SDS score.
tumor <- findOutliers(tumor,k=5,resolution=0.8)
tumor <- filter_cells_by_sds(tumor)
tumor <- tumor[,colData(tumor)$outlier== FALSE]

# KNN smoothing
tumor <- knnSmooth(tumor)

# scquantum & calcInteger
tumor <- calcInteger(tumor, method = 'scquantum', assay = 'smoothed_bincounts')

# clustering
tumor <- runUmap(tumor)
tumor <- findSuggestedK(tumor)
tumor <- findClusters(tumor)
tumor <- tumor[,colData(tumor)$subclones != 'c0']

# phylogenetic analysis
tumor <- runPhylo(tumor, metric = 'manhattan')

# calcConsensus
tumor <- calcConsensus(tumor)
```

## Visulization

```R
# add sample/visit information
colData(tumor)$Visit <- stringr::str_extract(colData(tumor)$sample, "(V0|V2)")
colData(tumor)$Info <- stringr::str_extract(colData(tumor)$sample, "\\d{4}_(V0|V2)")
colData(tumor)$Patient <- stringr::str_extract(colData(tumor)$sample, "\\d{4}")

# Umap
plotUmap(tumor, label = 'subclones')
plotUmap(tumor, label = 'Visit')

# Heatmap
# comparison between high/low quality cells
plotHeatmap(tumor, label = 'outlier', row_split ='outlier'

# consensus heatmap for each individual
# here you need to extract samples from the tumor first
patient <- tumor[,colData(tumor)$sample=='
plotHeatmap(tumor,consensus = TRUE,assay = 'integer',label = 'subclones',group = 'stage')

#

```
