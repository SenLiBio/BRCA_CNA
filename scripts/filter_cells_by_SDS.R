# Segment Deviation Score (SDS) and Correlation-Based Quality Control
# -------------------------------------------------------------------

# This script performs an additional layer of quality control on the CopyKit output.
# Specifically, it calculates the Segment Deviation Score (SDS), which quantifies the fit between the bin-level logR values and the corresponding segment-level ratios.

# For each bin, the squared difference between its logR ratio and the segment ratio it belongs to is calculated. 
# The average of these squared differences across all bins is reported as the SDS, which reflects the global smoothness and consistency of the profile.

# In addition, a KNN-based correlation score is computed for each sample. 
# Samples are retained only if SDS ≤ 1 and correlation ≥ 0.8 by default.
# These thresholds can be adjusted based on dataset-specific characteristics and logR profile inspection.

filter_cells_by_sds <- function(sample_obj) {
    segment_ratios_df <- as.data.frame(segment_ratios(sample_obj))
    ratios_df <- as.data.frame(ratios(sample_obj))
    if (!all(dim(segment_ratios_df) == dim(ratios_df))) {
        stop("Segment ratios and ratios must have the same dimensions.")
    }
    calculate_score <- function(segment_ratios, ratios) {
        diff_square <- (segment_ratios - ratios) ^ 2
        score <- sum(diff_square) / length(diff_square)
        return(score)
    }
    results <- sapply(1:ncol(segment_ratios_df), function(i) {
        calculate_score(segment_ratios_df[, i], ratios_df[, i])
    })
    results_df <- data.frame(cell = colnames(segment_ratios_df), score = results)
    results_df <- results_df[order(results_df$score), ]
    cell_names <- colnames(segment_ratios_df)
    if (!all(cell_names %in% rownames(colData(sample_obj)))) {
        stop("Mismatch between cell names and colData.")
    }
    colData(sample_obj)$score <- NA
    colData(sample_obj)$score[match(cell_names, rownames(colData(sample_obj)))] <- results
    unique_samples <- unique(colData(sample_obj)$sample)
    for (sample in unique_samples) {
        sample_indices <- which(colData(sample_obj)$sample == sample)
        colData(sample_obj)$outlier[sample_indices] <- ifelse(
            colData(sample_obj)$score[sample_indices] < 1 &
                colData(sample_obj)$cell_corr_value[sample_indices] >= 0.8,
            FALSE,
            TRUE
        )
    }
    sample_obj <- sample_obj[, colData(sample_obj)$outlier == FALSE]
    return(sample_obj)
}
