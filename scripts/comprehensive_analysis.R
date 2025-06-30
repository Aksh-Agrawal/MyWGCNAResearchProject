#!/usr/bin/env Rscript

# comprehensive_analysis.R
# ========================
# Complete WGCNA analysis pipeline with enhanced features

suppressPackageStartupMessages({
  library(optparse)
  library(MyWGCNAResearchProject)
  library(yaml)
  library(knitr)
  library(rmarkdown)
})

# Command line options
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/default.yaml",
              help = "Configuration file path"),
  make_option(c("-i", "--input"), type = "character", 
              help = "Input expression data file (overrides config)"),
  make_option(c("-o", "--output"), type = "character", default = "results",
              help = "Output directory"),
  make_option(c("--report"), action = "store_true", default = FALSE,
              help = "Generate HTML report"),
  make_option(c("--parallel"), action = "store_true", default = FALSE,
              help = "Enable parallel processing")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load configuration
config <- yaml::read_yaml(opt$config)

# Override input file if specified
if (!is.null(opt$input)) {
  config$expr_file <- opt$input
}

# Setup output directories
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "intermediate"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "reports"), recursive = TRUE, showWarnings = FALSE)

# Setup parallel processing
if (opt$parallel || config$use_parallel) {
  n_cores <- config$n_cores %||% parallel::detectCores() - 1
  WGCNA::enableWGCNAThreads(nThreads = n_cores)
  message(sprintf("Enabled parallel processing with %d cores", n_cores))
}

message("=== Starting Comprehensive WGCNA Analysis ===")

# Step 1: Load and preprocess data
message("\n1. Loading and preprocessing data...")

if (!file.exists(config$expr_file)) {
  # Create example data if input file doesn't exist
  message("Creating example expression data...")
  set.seed(123)
  n_genes <- 2000
  n_samples <- 30
  genes <- paste0("Gene_", 1:n_genes)
  samples <- paste0("Sample_", 1:n_samples)
  
  # Create expression with some modular structure
  expr_matrix <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 2), 
                       nrow = n_genes, ncol = n_samples,
                       dimnames = list(genes, samples))
  
  # Add modular structure
  module1_genes <- 1:200
  module2_genes <- 201:400
  module1_samples <- 1:15
  module2_samples <- 16:30
  
  expr_matrix[module1_genes, module1_samples] <- expr_matrix[module1_genes, module1_samples] + 3
  expr_matrix[module2_genes, module2_samples] <- expr_matrix[module2_genes, module2_samples] + 2.5
  
  expr_df <- data.frame(Gene = genes, expr_matrix, check.names = FALSE)
  dir.create(dirname(config$expr_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(expr_df, config$expr_file, row.names = FALSE)
  message(sprintf("Created example data: %s", config$expr_file))
}

# Load expression data
expr_raw <- loadExpressionData(config$expr_file)
message(sprintf("Loaded expression data: %d genes × %d samples", 
                nrow(expr_raw), ncol(expr_raw)))

# Advanced gene filtering
expr_filtered <- filterGenesByVariance(
  expr_raw,
  topN = config$topN,
  minExpr = config$min_expr,
  minSamples = config$min_samples,
  cv_threshold = config$cv_threshold
)

# Normalization
expr_norm <- normalizeExpression(expr_filtered, method = config$norm_method)
saveRDS(expr_norm, file.path(output_dir, "intermediate", "expr_normalized.rds"))

# Transpose for WGCNA (samples × genes)
datExpr <- t(expr_norm)

# Step 2: Advanced sample QC
message("\n2. Performing advanced sample quality control...")

qc_results <- advancedSampleQC(
  datExpr,
  zCut = config$outlier_z,
  mad_cut = config$mad_cut,
  connectivity_cut = config$connectivity_cut
)

datExpr_clean <- qc_results$cleanData
write.csv(qc_results$qcMetrics, 
         file.path(output_dir, "reports", "sample_qc_metrics.csv"),
         row.names = FALSE)

# Batch effect detection
if (config$detect_batch_effects) {
  batch_results <- detectBatchEffects(datExpr_clean)
  saveRDS(batch_results, file.path(output_dir, "intermediate", "batch_analysis.rds"))
}

# Step 3: Enhanced soft threshold selection
message("\n3. Enhanced soft threshold power selection...")

power_results <- pickSoftPower(
  datExpr_clean,
  powers = config$powers,
  networkType = config$tom_type,
  corType = config$cor_type,
  R2cut = config$R2_cut,
  slope_threshold = config$slope_threshold,
  mean_k_target = config$mean_k_target
)

selected_power <- power_results$power
saveRDS(power_results, file.path(output_dir, "intermediate", "power_selection.rds"))

# Plot power selection
png(file.path(output_dir, "figures", "soft_threshold_selection.png"),
    width = config$plot_width, height = config$plot_height, res = config$plot_dpi)
plotSoftPower(power_results$fitIndices, power_results$diagnostics, selected_power)
dev.off()

# Step 4: Advanced network construction
message("\n4. Constructing co-expression network...")

network_results <- constructNetwork(
  datExpr_clean,
  power = selected_power,
  networkType = config$tom_type,
  corType = config$cor_type,
  maxBlockSize = config$max_block_size,
  verbose = config$network_verbose
)

saveRDS(network_results$TOM, file.path(output_dir, "intermediate", "TOM_matrix.rds"))
saveRDS(network_results$network_metrics, 
        file.path(output_dir, "intermediate", "network_metrics.rds"))

# Step 5: Enhanced module detection
message("\n5. Detecting co-expression modules...")

module_results <- detectModules(
  network_results$dissTOM,
  deepSplit = config$deep_split,
  minModuleSize = config$min_module_size,
  method = "average",
  pamStage = config$pam_stage,
  alternativeMethods = config$alternative_methods
)

# Step 6: Module merging with validation
message("\n6. Merging similar modules...")

merge_results <- mergeModules(
  datExpr_clean,
  module_results$moduleColors,
  cutHeight = config$merge_height,
  minModuleSize = config$min_module_size,
  verbose = TRUE
)

final_colors <- merge_results$mergedColors
geneTree <- module_results$geneTree

# Save module assignments
module_df <- data.frame(
  Gene = colnames(datExpr_clean),
  Module = final_colors,
  stringsAsFactors = FALSE
)
write.table(module_df, 
           file.path(output_dir, "module_assignments.tsv"),
           sep = "\t", row.names = FALSE, quote = FALSE)

# Step 7: Enhanced visualization
message("\n7. Generating enhanced visualizations...")

# Dendrogram with modules
png(file.path(output_dir, "figures", "dendrogram_with_modules.png"),
    width = config$plot_width, height = config$plot_height, res = config$plot_dpi)
plotDendroAndColors(geneTree, final_colors,
                   groupLabels = "Modules",
                   main = "Gene Clustering Dendrogram with Module Colors",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
dev.off()

# Module eigengenes analysis
MEList <- moduleEigengenes(datExpr_clean, colors = final_colors)
MEs <- MEList$eigengenes

# Plot module eigengenes heatmap
if (ncol(MEs) > 1) {
  png(file.path(output_dir, "figures", "module_eigengenes_heatmap.png"),
      width = config$plot_width, height = config$plot_height, res = config$plot_dpi)
  ME_cor <- cor(MEs, use = "pairwise.complete.obs")
  heatmap(ME_cor, main = "Module Eigengene Correlation",
          col = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
}

# Step 8: Bootstrap stability assessment
if (config$bootstrap) {
  message("\n8. Assessing module stability via bootstrap...")
  
  stability_results <- replicate(config$n_boot, {
    # Bootstrap samples
    boot_idx <- sample(nrow(datExpr_clean), replace = TRUE)
    boot_expr <- datExpr_clean[boot_idx, ]
    
    # Construct bootstrap network
    boot_net <- constructNetwork(boot_expr, selected_power, 
                                config$tom_type, config$cor_type,
                                maxBlockSize = config$max_block_size,
                                verbose = FALSE)
    
    # Detect modules
    boot_modules <- detectModules(boot_net$dissTOM, 
                                 config$deep_split,
                                 config$min_module_size,
                                 pamStage = FALSE,
                                 alternativeMethods = c())
    
    # Compare with original modules
    # This is simplified - full implementation would use module preservation statistics
    length(unique(boot_modules$moduleColors)) - ("grey" %in% boot_modules$moduleColors)
  }, simplify = TRUE)
  
  saveRDS(stability_results, file.path(output_dir, "intermediate", "stability_assessment.rds"))
}

# Step 9: Generate comprehensive report
if (opt$report) {
  message("\n9. Generating comprehensive analysis report...")
  
  # Create report data
  report_data <- list(
    config = config,
    input_file = config$expr_file,
    n_samples = nrow(datExpr_clean),
    n_genes = ncol(datExpr_clean),
    selected_power = selected_power,
    power_diagnostics = power_results$diagnostics,
    n_modules = length(unique(final_colors)) - ("grey" %in% final_colors),
    module_sizes = table(final_colors)[names(table(final_colors)) != "grey"],
    network_metrics = network_results$network_metrics,
    qc_summary = qc_results$qcMetrics
  )
  
  saveRDS(report_data, file.path(output_dir, "intermediate", "report_data.rds"))
  
  # Generate simple text report for now
  report_file <- file.path(output_dir, "reports", "analysis_summary.txt")
  cat("WGCNA Analysis Summary\n", file = report_file)
  cat("======================\n\n", file = report_file, append = TRUE)
  cat(sprintf("Input: %s\n", config$expr_file), file = report_file, append = TRUE)
  cat(sprintf("Samples: %d\n", nrow(datExpr_clean)), file = report_file, append = TRUE)
  cat(sprintf("Genes: %d\n", ncol(datExpr_clean)), file = report_file, append = TRUE)
  cat(sprintf("Selected Power: %d\n", selected_power), file = report_file, append = TRUE)
  cat(sprintf("Modules Detected: %d\n", length(unique(final_colors)) - ("grey" %in% final_colors)), 
      file = report_file, append = TRUE)
  cat(sprintf("Unassigned Genes: %.1f%%\n", 
              100 * sum(final_colors == "grey") / length(final_colors)),
      file = report_file, append = TRUE)
}

message("\n=== Analysis Complete ===")
message(sprintf("Results saved to: %s", output_dir))
message(sprintf("Modules detected: %d", length(unique(final_colors)) - ("grey" %in% final_colors)))
message(sprintf("Selected power: %d", selected_power))

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x
