# soft_threshold.R
# =================
# Auto-pick and visualize soft-thresholding power for scale-free topology

library(WGCNA)

#' Advanced soft-thresholding power selection with multiple criteria
#'
#' @param datExpr samples × genes matrix
#' @param powers candidate integer vector
#' @param networkType "unsigned" or "signed"
#' @param corType "pearson" or "bicor"
#' @param R2cut minimum scale-free R²
#' @param slope_threshold minimum slope for scale-free fit
#' @param mean_k_target target mean connectivity
#' @return list(power, fitIndices, diagnostics)
pickSoftPower <- function(
    datExpr,
    powers = 1:20,
    networkType = "unsigned",
    corType = "pearson",
    R2cut = 0.80,
    slope_threshold = -1.0,
    mean_k_target = NULL
) {
  corFnc <- if (corType == "bicor") bicor else cor
  
  # Calculate soft threshold statistics
  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = networkType,
    corFnc = corFnc,
    verbose = 0
  )
  
  fit_indices <- sft$fitIndices
  
  # Enhanced power selection with multiple criteria
  
  # 1. Scale-free topology criterion
  R2_criterion <- fit_indices[, "SFT.R.sq"] >= R2cut
  
  # 2. Slope criterion (should be negative for scale-free)
  slope_criterion <- fit_indices[, "slope"] <= slope_threshold
  
  # 3. Mean connectivity criterion (avoid overly connected networks)
  if (is.null(mean_k_target)) {
    # Target connectivity ~ sqrt(number of samples)
    mean_k_target <- min(sqrt(nrow(datExpr)), 100)
  }
  mean_k <- fit_indices[, "mean.k."]
  connectivity_criterion <- mean_k <= mean_k_target & mean_k >= 1
  
  # 4. Combined criterion
  combined_criterion <- R2_criterion & slope_criterion & connectivity_criterion
  
  # Select power based on combined criteria
  power <- NA
  selection_reason <- "fallback"
  
  if (any(combined_criterion)) {
    # Choose minimum power that satisfies all criteria
    power <- min(powers[combined_criterion])
    selection_reason <- "combined_criteria"
  } else if (any(R2_criterion)) {
    # Fallback to R² criterion only
    power <- min(powers[R2_criterion])
    selection_reason <- "R2_only"
  } else {
    # Final fallback to highest R²
    power <- powers[which.max(fit_indices[, "SFT.R.sq"])]
    selection_reason <- "max_R2"
  }
  
  # Diagnostics
  selected_idx <- which(powers == power)
  diagnostics <- list(
    selected_power = power,
    selected_R2 = fit_indices[selected_idx, "SFT.R.sq"],
    selected_slope = fit_indices[selected_idx, "slope"],
    selected_mean_k = fit_indices[selected_idx, "mean.k."],
    selection_reason = selection_reason,
    powers_meeting_R2 = powers[R2_criterion],
    powers_meeting_slope = powers[slope_criterion],
    powers_meeting_connectivity = powers[connectivity_criterion],
    target_mean_k = mean_k_target
  )
  
  message(sprintf("Selected power %d (R²=%.3f, slope=%.2f, mean.k=%.1f) - %s",
                  power, diagnostics$selected_R2, diagnostics$selected_slope,
                  diagnostics$selected_mean_k, selection_reason))
  
  list(power = power, fitIndices = fit_indices, diagnostics = diagnostics)
}

#' Enhanced visualization of soft power selection with diagnostic plots
#'
#' @param fitIndices from pickSoftThreshold()
#' @param diagnostics diagnostic information from enhanced power selection
#' @param selected_power the chosen power
plotSoftPower <- function(fitIndices, diagnostics = NULL, selected_power = NULL) {
  
  # Set up multi-panel plot
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  powers <- fitIndices[, 1]
  
  # 1. Scale-free topology fit
  plot(powers, fitIndices[, 2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R²",
       type = "b", main = "Scale-free topology fit")
  abline(h = 0.8, col = "red", lty = 2)
  if (!is.null(selected_power)) {
    abline(v = selected_power, col = "blue", lty = 2)
    points(selected_power, fitIndices[powers == selected_power, 2], 
           col = "blue", pch = 19, cex = 1.5)
  }
  
  # 2. Mean connectivity
  plot(powers, fitIndices[, 5],
       xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity",
       type = "b", main = "Mean connectivity")
  if (!is.null(diagnostics) && !is.null(diagnostics$target_mean_k)) {
    abline(h = diagnostics$target_mean_k, col = "red", lty = 2)
  }
  if (!is.null(selected_power)) {
    abline(v = selected_power, col = "blue", lty = 2)
    points(selected_power, fitIndices[powers == selected_power, 5], 
           col = "blue", pch = 19, cex = 1.5)
  }
  
  # 3. Slope of scale-free fit
  if (ncol(fitIndices) >= 6) {
    plot(powers, fitIndices[, 6],
         xlab = "Soft Threshold (power)", 
         ylab = "Slope",
         type = "b", main = "Scale-free fit slope")
    abline(h = -1, col = "red", lty = 2)
    if (!is.null(selected_power)) {
      abline(v = selected_power, col = "blue", lty = 2)
      points(selected_power, fitIndices[powers == selected_power, 6], 
             col = "blue", pch = 19, cex = 1.5)
    }
  }
  
  # 4. Summary statistics
  plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "", ylab = "", main = "Selection Summary", axes = FALSE)
  
  if (!is.null(diagnostics)) {
    summary_text <- paste(
      sprintf("Selected Power: %d", diagnostics$selected_power),
      sprintf("R²: %.3f", diagnostics$selected_R2),
      sprintf("Slope: %.2f", diagnostics$selected_slope),
      sprintf("Mean k: %.1f", diagnostics$selected_mean_k),
      sprintf("Method: %s", diagnostics$selection_reason),
      sep = "\n"
    )
    text(0.5, 0.5, summary_text, cex = 0.9, adj = 0.5)
  }
  
  par(mfrow = c(1, 1))  # Reset layout
}
