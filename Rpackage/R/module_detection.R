# module_detection.R
# ====================
# Dynamic tree cutting & module merging

library(WGCNA)
library(dynamicTreeCut)

#' Advanced module detection with multiple algorithms and validation
#'
#' @param dissTOM TOM dissimilarity matrix
#' @param deepSplit sensitivity (0–4)
#' @param minModuleSize minimum cluster size
#' @param method clustering linkage ("average", etc.)
#' @param pamStage whether to use PAM stage in dynamic tree cutting
#' @param alternativeMethods additional clustering methods to compare
#' @return list(geneTree, moduleColors, methods_comparison, stability_metrics)
detectModules <- function(
    dissTOM,
    deepSplit = 2,
    minModuleSize = 30,
    method = "average",
    pamStage = TRUE,
    alternativeMethods = c("ward.D2", "complete")
) {
  
  # Primary clustering with specified method
  geneTree <- flashClust(as.dist(dissTOM), method = method)
  
  # Dynamic tree cutting with enhanced parameters
  dynMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = deepSplit,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize,
    pamStage = pamStage,
    verbose = 0
  )
  moduleColors <- labels2colors(dynMods)
  
  # Compare with alternative clustering methods
  methods_comparison <- list()
  
  for (alt_method in alternativeMethods) {
    tryCatch({
      alt_tree <- flashClust(as.dist(dissTOM), method = alt_method)
      alt_mods <- cutreeDynamic(
        dendro = alt_tree,
        distM = dissTOM,
        deepSplit = deepSplit,
        pamRespectsDendro = FALSE,
        minClusterSize = minModuleSize,
        verbose = 0
      )
      alt_colors <- labels2colors(alt_mods)
      
      # Calculate agreement with primary method
      agreement <- sum(dynMods == alt_mods) / length(dynMods)
      
      methods_comparison[[alt_method]] <- list(
        tree = alt_tree,
        modules = alt_mods,
        colors = alt_colors,
        agreement_with_primary = agreement,
        num_modules = length(unique(alt_colors)) - ("grey" %in% alt_colors)
      )
    }, error = function(e) {
      warning(sprintf("Alternative method %s failed: %s", alt_method, e$message))
    })
  }
  
  # Calculate stability metrics
  stability_metrics <- list(
    num_modules = length(unique(moduleColors)) - ("grey" %in% moduleColors),
    largest_module_size = max(table(moduleColors)),
    grey_fraction = sum(moduleColors == "grey") / length(moduleColors),
    cophenetic_correlation = cor(cophenetic(geneTree), as.dist(dissTOM), use = "complete.obs")
  )
  
  # Module size distribution
  module_sizes <- table(moduleColors)
  module_sizes <- module_sizes[names(module_sizes) != "grey"]
  stability_metrics$module_size_stats <- list(
    mean_size = mean(module_sizes),
    median_size = median(module_sizes),
    size_cv = sd(module_sizes) / mean(module_sizes)
  )
  
  message(sprintf("Detected %d modules (%.1f%% unassigned genes)", 
                  stability_metrics$num_modules, 
                  100 * stability_metrics$grey_fraction))
  
  list(
    geneTree = geneTree, 
    moduleColors = moduleColors,
    methods_comparison = methods_comparison,
    stability_metrics = stability_metrics
  )
}

#' Enhanced module merging with multiple criteria and validation
#'
#' @param datExpr samples × genes matrix
#' @param moduleColors vector of module assignments
#' @param cutHeight merge threshold (e.g. 0.25)
#' @param minModuleSize minimum size for modules after merging
#' @param verbose whether to print merging details
#' @return list(mergedColors, mergeInfo, merge_statistics)
mergeModules <- function(
    datExpr,
    moduleColors,
    cutHeight = 0.25,
    minModuleSize = 30,
    verbose = TRUE
) {
  
  # Store original module information
  original_modules <- unique(moduleColors)
  original_count <- length(original_modules) - ("grey" %in% original_modules)
  
  # Enhanced module merging
  merge <- mergeCloseModules(
    datExpr,
    moduleColors,
    cutHeight = cutHeight,
    verbose = if(verbose) 2 else 0
  )
  
  mergedColors <- merge$colors
  
  # Post-merge cleanup: reassign small modules to grey
  module_sizes <- table(mergedColors)
  small_modules <- names(module_sizes)[module_sizes < minModuleSize & names(module_sizes) != "grey"]
  
  if (length(small_modules) > 0) {
    mergedColors[mergedColors %in% small_modules] <- "grey"
    if (verbose) {
      message(sprintf("Reassigned %d small modules to grey: %s", 
                      length(small_modules), 
                      paste(small_modules, collapse = ", ")))
    }
  }
  
  # Calculate merge statistics
  final_modules <- unique(mergedColors)
  final_count <- length(final_modules) - ("grey" %in% final_modules)
  
  merge_statistics <- list(
    original_module_count = original_count,
    final_module_count = final_count,
    modules_merged = original_count - final_count,
    merge_efficiency = (original_count - final_count) / original_count,
    final_grey_fraction = sum(mergedColors == "grey") / length(mergedColors),
    size_distribution = table(mergedColors)[names(table(mergedColors)) != "grey"]
  )
  
  # Validate merged modules
  validation_metrics <- list()
  
  if (final_count > 0) {
    # Calculate module eigengenes for validation
    MEList <- moduleEigengenes(datExpr, colors = mergedColors)
    MEs <- MEList$eigengenes
    
    # Module preservation statistics
    if (ncol(MEs) > 1) {
      ME_cor <- cor(MEs, use = "pairwise.complete.obs")
      diag(ME_cor) <- NA
      validation_metrics$max_inter_module_correlation <- max(abs(ME_cor), na.rm = TRUE)
      validation_metrics$mean_inter_module_correlation <- mean(abs(ME_cor), na.rm = TRUE)
    }
    
    # Module coherence (average intra-module correlation)
    module_coherence <- sapply(final_modules[final_modules != "grey"], function(mod) {
      mod_genes <- names(mergedColors)[mergedColors == mod]
      if (length(mod_genes) >= 2) {
        mod_expr <- datExpr[, mod_genes, drop = FALSE]
        mod_cor <- cor(mod_expr, use = "pairwise.complete.obs")
        mean(mod_cor[upper.tri(mod_cor)], na.rm = TRUE)
      } else {
        NA
      }
    })
    
    validation_metrics$module_coherence <- module_coherence
    validation_metrics$mean_module_coherence <- mean(module_coherence, na.rm = TRUE)
  }
  
  if (verbose) {
    message(sprintf("Module merging: %d -> %d modules (%.1f%% reduction)", 
                    original_count, final_count, 
                    100 * merge_statistics$merge_efficiency))
  }
  
  list(
    mergedColors = mergedColors, 
    mergeInfo = merge,
    merge_statistics = merge_statistics,
    validation_metrics = validation_metrics
  )
}
