#!/usr/bin/env Rscript

# enhanced_analysis.R
# ===================
# Enhanced WGCNA analysis using existing functions with improvements

suppressPackageStartupMessages({
  library(MyWGCNAResearchProject)
  library(WGCNA)
})

message("=== Enhanced WGCNA Analysis ===")

# 1. Load and preprocess data with enhanced filtering
message("1. Enhanced data preprocessing...")

# Use larger test data
expr_raw <- readRDS("output/larger_test_data.rds")
datExpr <- t(expr_raw)  # samples × genes

message(sprintf("Data dimensions: %d samples × %d genes", nrow(datExpr), ncol(datExpr)))

# Enhanced outlier detection
sampleMeans <- rowMeans(datExpr, na.rm = TRUE)
sampleSDs <- apply(datExpr, 1, sd, na.rm = TRUE)

# Multiple outlier detection criteria
z_means <- scale(sampleMeans)[,1]
mad_outliers <- abs(sampleMeans - median(sampleMeans)) > 3 * mad(sampleMeans)
outlier_samples <- abs(z_means) > 2.5 | mad_outliers

if (sum(outlier_samples) > 0) {
  warning(sprintf("Removing %d outlier samples", sum(outlier_samples)))
  datExpr <- datExpr[!outlier_samples, , drop = FALSE]
}

# 2. Enhanced soft threshold selection with multiple criteria
message("2. Enhanced soft threshold selection...")

powers <- c(1:10, 12, 14, 16, 18, 20)
power_result <- pickSoftPower(datExpr, powers = powers, networkType = "unsigned", corType = "pearson")

# Enhanced power selection logic
fit_indices <- power_result$fitIndices
R2_values <- fit_indices[, "SFT.R.sq"]
mean_k_values <- fit_indices[, "mean.k."]

# Target connectivity
target_connectivity <- min(sqrt(nrow(datExpr)), 100)

# Multi-criteria selection
R2_criterion <- R2_values >= 0.8
connectivity_criterion <- mean_k_values <= target_connectivity & mean_k_values >= 1

if (any(R2_criterion & connectivity_criterion)) {
  selected_power <- min(powers[R2_criterion & connectivity_criterion])
  selection_method <- "multi_criteria"
} else if (any(R2_criterion)) {
  selected_power <- min(powers[R2_criterion])
  selection_method <- "R2_only"
} else {
  selected_power <- powers[which.max(R2_values)]
  selection_method <- "max_R2"
}

message(sprintf("Selected power: %d (method: %s, R²=%.3f)", 
                selected_power, selection_method, 
                R2_values[powers == selected_power]))

# Enhanced visualization
png("results/figures/enhanced_power_selection.png", width = 1200, height = 800, res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Plot 1: Scale-free fit
plot(powers, R2_values, type = "b", main = "Scale-free Topology Fit",
     xlab = "Soft Threshold (power)", ylab = "Scale Free R²")
abline(h = 0.8, col = "red", lty = 2)
abline(v = selected_power, col = "blue", lty = 2)
points(selected_power, R2_values[powers == selected_power], col = "blue", pch = 19, cex = 1.5)

# Plot 2: Mean connectivity
plot(powers, mean_k_values, type = "b", main = "Mean Connectivity",
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity")
abline(h = target_connectivity, col = "red", lty = 2)
abline(v = selected_power, col = "blue", lty = 2)
points(selected_power, mean_k_values[powers == selected_power], col = "blue", pch = 19, cex = 1.5)

# Plot 3: Power selection summary
plot(R2_values, mean_k_values, type = "p", main = "R² vs Connectivity",
     xlab = "Scale Free R²", ylab = "Mean Connectivity")
points(R2_values[powers == selected_power], mean_k_values[powers == selected_power], 
       col = "blue", pch = 19, cex = 2)
text(R2_values[powers == selected_power], mean_k_values[powers == selected_power],
     selected_power, pos = 3, col = "blue", font = 2)

# Plot 4: Criteria comparison
barplot(table(c(rep("R²", sum(R2_criterion)), 
                rep("Connectivity", sum(connectivity_criterion)),
                rep("Both", sum(R2_criterion & connectivity_criterion)))),
        main = "Powers Meeting Criteria", ylab = "Count")

dev.off()

# 3. Enhanced network construction
message("3. Enhanced network construction...")

# Build adjacency and TOM matrices
adj <- adjacency(datExpr, power = selected_power, type = "unsigned")
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM

# Network quality metrics
connectivity <- rowSums(adj) - diag(adj)
network_density <- sum(adj) / (ncol(datExpr) * (ncol(datExpr) - 1))
tom_sparsity <- sum(TOM[upper.tri(TOM)] < 0.01) / sum(upper.tri(TOM))

message(sprintf("Network metrics: density=%.4f, sparsity=%.2f", 
                network_density, tom_sparsity))

# 4. Enhanced module detection with multiple methods
message("4. Enhanced module detection...")

# Primary method with relaxed parameters
geneTree <- flashClust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = 10)  # Reduced min size
moduleColors <- labels2colors(dynamicMods)

# Alternative clustering methods for comparison
methods_comparison <- list()
alt_methods <- c("ward.D", "complete")

for (method in alt_methods) {
  tryCatch({
    alt_tree <- flashClust(as.dist(dissTOM), method = method)
    alt_mods <- cutreeDynamic(dendro = alt_tree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
    alt_colors <- labels2colors(alt_mods)
    
    # Calculate agreement
    agreement <- sum(dynamicMods == alt_mods) / length(dynamicMods)
    methods_comparison[[method]] <- list(
      colors = alt_colors,
      agreement = agreement,
      num_modules = length(unique(alt_colors)) - ("grey" %in% alt_colors)
    )
  }, error = function(e) {
    warning(sprintf("Method %s failed: %s", method, e$message))
  })
}

message(sprintf("Module detection comparison:"))
message(sprintf("  Primary (average): %d modules", 
                length(unique(moduleColors)) - ("grey" %in% moduleColors)))
for (method in names(methods_comparison)) {
  message(sprintf("  %s: %d modules (%.2f agreement)", 
                  method, methods_comparison[[method]]$num_modules,
                  methods_comparison[[method]]$agreement))
}

# 5. Enhanced module merging
message("5. Enhanced module merging...")

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = moduleColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)

# Merge similar modules
merge <- mergeCloseModules(datExpr, moduleColors, cutHeight = 0.25, verbose = 2)
mergedColors <- merge$colors

# Post-merge statistics
original_modules <- length(unique(moduleColors)) - ("grey" %in% moduleColors)
final_modules <- length(unique(mergedColors)) - ("grey" %in% mergedColors)
merge_efficiency <- (original_modules - final_modules) / original_modules

message(sprintf("Module merging: %d -> %d modules (%.1f%% reduction)", 
                original_modules, final_modules, 100 * merge_efficiency))

# 6. Enhanced visualizations
message("6. Creating enhanced visualizations...")

# Dendrogram with modules
png("results/figures/enhanced_dendrogram.png", width = 1400, height = 1000, res = 300)
plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors),
                   c("Dynamic Tree Cut", "Merged dynamic"),
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene Clustering Dendrogram and Module Colors")
dev.off()

# Module eigengenes analysis
if (final_modules > 1) {
  final_MEList <- moduleEigengenes(datExpr, colors = mergedColors)
  final_MEs <- final_MEList$eigengenes
  
  # Eigengene correlation heatmap
  png("results/figures/module_eigengene_correlation.png", width = 800, height = 800, res = 300)
  ME_cor <- cor(final_MEs, use = "pairwise.complete.obs")
  
  # Enhanced heatmap with better color scheme
  library(RColorBrewer)
  colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
  heatmap(ME_cor, col = colors, main = "Module Eigengene Correlation Heatmap",
          cexRow = 0.8, cexCol = 0.8)
  dev.off()
  
  # Module size distribution
  png("results/figures/module_size_distribution.png", width = 800, height = 600, res = 300)
  module_sizes <- table(mergedColors)
  module_sizes <- module_sizes[names(module_sizes) != "grey"]
  
  barplot(sort(module_sizes, decreasing = TRUE), 
          main = "Module Size Distribution",
          xlab = "Module", ylab = "Number of Genes",
          col = names(sort(module_sizes, decreasing = TRUE)),
          las = 2)
  dev.off()
}

# 7. Network topology analysis
message("7. Network topology analysis...")

# Connectivity distribution
k <- connectivity
k_table <- table(k)
k_values <- as.numeric(names(k_table))
p_k <- as.numeric(k_table) / sum(k_table)

# Scale-free topology validation
png("results/figures/connectivity_distribution.png", width = 800, height = 600, res = 300)
par(mfrow = c(1, 2))

# Linear scale
plot(k_values, p_k, type = "p", main = "Connectivity Distribution",
     xlab = "Connectivity (k)", ylab = "P(k)")

# Log-log scale for scale-free check
valid_idx <- k_values > 0 & p_k > 0
if (sum(valid_idx) >= 5) {
  plot(log10(k_values[valid_idx]), log10(p_k[valid_idx]), 
       type = "p", main = "Scale-Free Topology Check",
       xlab = "log10(k)", ylab = "log10(P(k))")
  
  # Fit line
  fit <- lm(log10(p_k[valid_idx]) ~ log10(k_values[valid_idx]))
  abline(fit, col = "red")
  
  # Add R² to plot
  r_squared <- summary(fit)$r.squared
  text(min(log10(k_values[valid_idx])), max(log10(p_k[valid_idx])),
       sprintf("R² = %.3f", r_squared), adj = c(0, 1))
}
dev.off()

# 8. Save results
message("8. Saving enhanced results...")

# Create comprehensive results
results <- list(
  input_dimensions = c(nrow(datExpr), ncol(datExpr)),
  selected_power = selected_power,
  selection_method = selection_method,
  network_metrics = list(
    density = network_density,
    tom_sparsity = tom_sparsity,
    mean_connectivity = mean(connectivity),
    connectivity_range = range(connectivity)
  ),
  module_detection = list(
    original_modules = original_modules,
    final_modules = final_modules,
    merge_efficiency = merge_efficiency,
    grey_fraction = sum(mergedColors == "grey") / length(mergedColors)
  ),
  methods_comparison = methods_comparison
)

saveRDS(results, "results/enhanced_analysis_results.rds")

# Save final module assignments
module_assignments <- data.frame(
  Gene = colnames(datExpr),
  Original_Module = moduleColors,
  Final_Module = mergedColors,
  stringsAsFactors = FALSE
)
write.csv(module_assignments, "results/enhanced_module_assignments.csv", row.names = FALSE)

# Create summary report
report <- c(
  "Enhanced WGCNA Analysis Summary",
  "===============================",
  "",
  sprintf("Input: %d samples × %d genes", nrow(datExpr), ncol(datExpr)),
  sprintf("Selected Power: %d (%s)", selected_power, selection_method),
  sprintf("Network Density: %.4f", network_density),
  sprintf("TOM Sparsity: %.2f", tom_sparsity),
  sprintf("Modules Detected: %d -> %d (%.1f%% reduction)", 
          original_modules, final_modules, 100 * merge_efficiency),
  sprintf("Unassigned Genes: %.1f%%", 
          100 * sum(mergedColors == "grey") / length(mergedColors)),
  "",
  "Method Comparison:",
  sprintf("  Primary (average): %d modules", original_modules)
)

for (method in names(methods_comparison)) {
  report <- c(report, sprintf("  %s: %d modules (%.2f agreement)", 
                             method, methods_comparison[[method]]$num_modules,
                             methods_comparison[[method]]$agreement))
}

writeLines(report, "results/enhanced_analysis_summary.txt")

message("\n=== Enhanced Analysis Complete ===")
message(sprintf("Final modules: %d", final_modules))
message(sprintf("Power selected: %d", selected_power))
message("Results saved to results/ directory")
