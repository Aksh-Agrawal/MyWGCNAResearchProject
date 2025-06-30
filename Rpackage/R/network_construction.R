# network_construction.R
# =======================
# Build adjacency and TOM matrices

library(WGCNA)

#' Advanced network construction with multiple validation metrics
#'
#' @param datExpr samples × genes matrix
#' @param power soft-thresholding power
#' @param networkType "unsigned" or "signed"
#' @param corType "pearson" or "bicor"
#' @param maxBlockSize maximum genes per block for large datasets
#' @param verbose whether to print progress
#' @return list(adjacency, TOM, dissTOM, network_metrics, block_info)
constructNetwork <- function(
    datExpr,
    power,
    networkType = "unsigned",
    corType = "pearson",
    maxBlockSize = 5000,
    verbose = TRUE
) {
  
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  if (verbose) {
    message(sprintf("Constructing network: %d genes, %d samples, power = %d", 
                    nGenes, nSamples, power))
  }
  
  # Select correlation function
  corFnc <- if (corType == "bicor") bicor else cor
  
  # Handle large datasets with block processing
  if (nGenes > maxBlockSize) {
    if (verbose) message("Using block-wise network construction for large dataset")
    
    # Block-wise network construction
    net <- blockwiseModules(
      datExpr,
      power = power,
      networkType = networkType,
      corType = corType,
      maxBlockSize = maxBlockSize,
      TOMType = networkType,
      saveTOMs = FALSE,
      verbose = if(verbose) 2 else 0,
      returnTOMs = TRUE
    )
    
    # Extract components
    TOM <- net$TOMFiles[[1]]  # For single block or combined
    adjacency <- NULL  # Block-wise doesn't return full adjacency
    dissTOM <- 1 - TOM
    
    block_info <- list(
      num_blocks = length(net$blocks),
      block_sizes = sapply(net$blocks, length),
      block_method = "blockwise"
    )
    
  } else {
    # Standard single-block construction
    if (verbose) message("Using standard network construction")
    
    # Calculate adjacency matrix
    adjacency <- adjacency(
      datExpr,
      power = power,
      type = networkType,
      corFnc = corFnc
    )
    
    # Calculate TOM matrix
    TOM <- TOMsimilarity(adjacency, TOMType = networkType)
    dissTOM <- 1 - TOM
    
    block_info <- list(
      num_blocks = 1,
      block_sizes = nGenes,
      block_method = "single"
    )
  }
  
  # Calculate network metrics
  network_metrics <- list()
  
  # Basic network properties
  if (!is.null(adjacency)) {
    # Connectivity statistics
    connectivity <- rowSums(adjacency) - diag(adjacency)
    network_metrics$connectivity_stats <- list(
      mean = mean(connectivity),
      median = median(connectivity),
      sd = sd(connectivity),
      range = range(connectivity)
    )
    
    # Network density
    network_metrics$density <- sum(adjacency) / (nGenes * (nGenes - 1))
    
    # Clustering coefficient (approximate for large networks)
    if (nGenes <= 1000) {
      clustering_coef <- sapply(1:nGenes, function(i) {
        neighbors <- which(adjacency[i, ] > 0.1)  # threshold for "connected"
        if (length(neighbors) < 2) return(0)
        subgraph <- adjacency[neighbors, neighbors]
        sum(subgraph > 0.1) / (length(neighbors) * (length(neighbors) - 1))
      })
      network_metrics$mean_clustering_coefficient <- mean(clustering_coef, na.rm = TRUE)
    }
  }
  
  # TOM-based metrics
  if (!is.null(TOM)) {
    tom_values <- TOM[upper.tri(TOM)]
    network_metrics$TOM_stats <- list(
      mean = mean(tom_values, na.rm = TRUE),
      median = median(tom_values, na.rm = TRUE),
      quantiles = quantile(tom_values, c(0.05, 0.25, 0.75, 0.95), na.rm = TRUE)
    )
    
    # Sparsity
    network_metrics$TOM_sparsity <- sum(tom_values < 0.01) / length(tom_values)
  }
  
  # Scale-free topology validation
  if (!is.null(adjacency) && nGenes <= 2000) {  # Only for manageable sizes
    k <- rowSums(adjacency > 0.1)  # degree with threshold
    if (length(unique(k)) > 5) {  # Need enough different degrees
      k_table <- table(k)
      k_values <- as.numeric(names(k_table))
      p_k <- as.numeric(k_table) / sum(k_table)
      
      # Fit power law: log(p(k)) ~ -gamma * log(k)
      valid_idx <- k_values > 0 & p_k > 0
      if (sum(valid_idx) >= 5) {
        log_k <- log10(k_values[valid_idx])
        log_p <- log10(p_k[valid_idx])
        
        fit <- lm(log_p ~ log_k)
        network_metrics$scale_free_fit <- list(
          gamma = -coef(fit)[2],  # Power law exponent
          R_squared = summary(fit)$r.squared,
          p_value = summary(fit)$coefficients[2, 4]
        )
      }
    }
  }
  
  if (verbose) {
    message(sprintf("Network construction complete. TOM sparsity: %.2f", 
                    network_metrics$TOM_sparsity %||% NA))
  }
  
  list(
    adjacency = adjacency, 
    TOM = TOM, 
    dissTOM = dissTOM,
    network_metrics = network_metrics,
    block_info = block_info
  )
}

# Helper operator for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x
