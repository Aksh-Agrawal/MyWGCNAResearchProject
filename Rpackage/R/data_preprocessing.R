# data_preprocessing.R
# ====================
# Functions to load, QC-filter and normalize expression data

# dependencies
library(WGCNA)
if (!requireNamespace("DESeq2", quietly=TRUE)) {
  message("installing DESeq2 for VST normalization...")
  BiocManager::install("DESeq2")
}
library(DESeq2)
library(S4Vectors)

#' Load expression matrix (genes × samples)
#'
#' @param file path to CSV/TSV
#' @param sep field separator (',' or '\t')
#' @param row.names column for gene IDs
#' @return numeric matrix
loadExpressionData <- function(file, sep = ",", row.names = 1) {
  df <- read.csv(file, sep = sep, row.names = row.names, check.names = FALSE)
  as.matrix(df)
}

#' Advanced gene filtering with multiple criteria
#'
#' @param exprData genes × samples matrix
#' @param topN number of highest-variance genes to keep
#' @param minExpr minimum expression threshold
#' @param minSamples minimum samples with expression above threshold
#' @param cv_threshold coefficient of variation threshold
#' @return filtered matrix with filtering statistics
filterGenesByVariance <- function(exprData, topN = 10000, minExpr = 1, 
                                  minSamples = 3, cv_threshold = 0.1) {
  # Record initial dimensions
  initial_genes <- nrow(exprData)
  
  # 1. Remove low-expressed genes
  expressed <- rowSums(exprData >= minExpr) >= minSamples
  exprData <- exprData[expressed, , drop = FALSE]
  
  # 2. Remove genes with very low variability
  vars <- apply(exprData, 1, var, na.rm = TRUE)
  means <- rowMeans(exprData, na.rm = TRUE)
  cv <- sqrt(vars) / means
  cv[is.infinite(cv) | is.nan(cv)] <- 0
  
  # Keep genes with CV above threshold
  variable_genes <- cv >= cv_threshold
  exprData <- exprData[variable_genes, , drop = FALSE]
  
  # 3. Select top N most variable genes
  vars_filtered <- vars[variable_genes]
  if (nrow(exprData) > topN) {
    keep <- order(vars_filtered, decreasing = TRUE)[seq_len(topN)]
    exprData <- exprData[keep, , drop = FALSE]
  }
  
  # Print filtering statistics
  message(sprintf("Gene filtering: %d -> %d genes (%.1f%% retained)", 
                  initial_genes, nrow(exprData), 
                  100 * nrow(exprData) / initial_genes))
  
  attr(exprData, "filtering_stats") <- list(
    initial_genes = initial_genes,
    final_genes = nrow(exprData),
    expressed_genes = sum(expressed),
    variable_genes = sum(variable_genes)
  )
  
  exprData
}

#' Detect & remove outlier samples by mean-expression Z-score
#'
#' @param datExpr samples × genes matrix
#' @param zCut Z-score cutoff
#' @return cleaned datExpr
detectOutlierSamples <- function(datExpr, zCut = 2.5) {
  sampleMeans <- rowMeans(datExpr, na.rm = TRUE)
  z <- scale(sampleMeans)
  out <- which(abs(z) > zCut)
  if (length(out)) {
    warning("Removing outlier samples: ", paste(rownames(datExpr)[out], collapse = ", "))
    datExpr <- datExpr[-out, , drop = FALSE]
  }
  datExpr
}

#' Advanced sample quality control with multiple metrics
#'
#' @param datExpr samples × genes matrix
#' @param zCut Z-score cutoff for outlier detection
#' @param mad_cut MAD-based outlier detection threshold
#' @param connectivity_cut Network connectivity threshold
#' @return list with cleaned data and QC metrics
advancedSampleQC <- function(datExpr, zCut = 2.5, mad_cut = 3, connectivity_cut = -2.5) {
  
  # Calculate multiple QC metrics
  sampleMeans <- rowMeans(datExpr, na.rm = TRUE)
  sampleSDs <- apply(datExpr, 1, sd, na.rm = TRUE)
  sampleMedians <- apply(datExpr, 1, median, na.rm = TRUE)
  
  # 1. Mean-based outlier detection
  z_means <- scale(sampleMeans)[,1]
  outliers_mean <- abs(z_means) > zCut
  
  # 2. MAD-based outlier detection  
  mad_val <- mad(sampleMeans)
  median_val <- median(sampleMeans)
  mad_outliers <- abs(sampleMeans - median_val) > mad_cut * mad_val
  
  # 3. Network connectivity-based outlier detection
  if (ncol(datExpr) >= 50) {  # Only for sufficient genes
    sampleTree <- hclust(dist(datExpr), method = "average")
    connectivity <- fundamentalNetworkConcepts(adjacency(datExpr))$Connectivity
    z_connectivity <- scale(connectivity)[,1]
    connectivity_outliers <- z_connectivity < connectivity_cut
  } else {
    connectivity_outliers <- rep(FALSE, nrow(datExpr))
    z_connectivity <- rep(0, nrow(datExpr))
  }
  
  # 4. Combine outlier detection methods
  all_outliers <- outliers_mean | mad_outliers | connectivity_outliers
  outlier_samples <- rownames(datExpr)[all_outliers]
  
  # Remove outliers
  if (sum(all_outliers) > 0) {
    warning(sprintf("Removing %d outlier samples: %s", 
                    sum(all_outliers), 
                    paste(outlier_samples, collapse = ", ")))
    datExpr_clean <- datExpr[!all_outliers, , drop = FALSE]
  } else {
    datExpr_clean <- datExpr
  }
  
  # QC metrics
  qc_metrics <- data.frame(
    Sample = rownames(datExpr),
    Mean = sampleMeans,
    SD = sampleSDs,
    Median = sampleMedians,
    Z_Mean = z_means,
    MAD_Outlier = mad_outliers,
    Connectivity = if(length(z_connectivity) == nrow(datExpr)) z_connectivity else NA,
    Is_Outlier = all_outliers,
    stringsAsFactors = FALSE
  )
  
  list(
    cleanData = datExpr_clean,
    qcMetrics = qc_metrics,
    removedSamples = outlier_samples
  )
}

#' Detect and visualize batch effects
#'
#' @param datExpr samples × genes matrix
#' @param batch_info data.frame with sample metadata including batch
#' @return list with PCA results and batch effect statistics
detectBatchEffects <- function(datExpr, batch_info = NULL) {
  
  # PCA analysis
  pca_result <- prcomp(datExpr, center = TRUE, scale. = TRUE)
  variance_explained <- summary(pca_result)$importance[2,]
  
  # If batch information available, test for batch effects
  batch_stats <- NULL
  if (!is.null(batch_info) && "batch" %in% colnames(batch_info)) {
    # Match samples
    common_samples <- intersect(rownames(datExpr), rownames(batch_info))
    if (length(common_samples) > 0) {
      pc1_scores <- pca_result$x[common_samples, "PC1"]
      pc2_scores <- pca_result$x[common_samples, "PC2"]
      batch_labels <- batch_info[common_samples, "batch"]
      
      # ANOVA test for batch effects
      batch_p_pc1 <- tryCatch({
        anova(lm(pc1_scores ~ as.factor(batch_labels)))$`Pr(>F)`[1]
      }, error = function(e) NA)
      
      batch_p_pc2 <- tryCatch({
        anova(lm(pc2_scores ~ as.factor(batch_labels)))$`Pr(>F)`[1]
      }, error = function(e) NA)
      
      batch_stats <- list(
        pc1_batch_p = batch_p_pc1,
        pc2_batch_p = batch_p_pc2,
        significant_batch_effect = any(c(batch_p_pc1, batch_p_pc2) < 0.05, na.rm = TRUE)
      )
    }
  }
  
  list(
    pca = pca_result,
    variance_explained = variance_explained,
    batch_stats = batch_stats
  )
}

#' Normalize using log2 or variance-stabilizing transform (VST)
#'
#' @param exprData raw counts genes × samples
#' @param method "log2" or "vst"
#' @return normalized matrix
normalizeExpression <- function(exprData, method = c("log2", "vst")) {
  method <- match.arg(method)
  if (method == "log2") {
    log2(exprData + 1)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = exprData,
                                  colData = DataFrame(row = colnames(exprData)),
                                  design = ~ 1)
    vst(dds, blind = TRUE) |> assay()
  }
}
