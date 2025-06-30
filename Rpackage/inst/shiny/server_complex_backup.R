# Premium Enhanced Server for WGCNA Explorer
# =============================================

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(MyWGCNAResearchProject)
library(plotly)
library(DT)
library(WGCNA)
library(RColorBrewer)
library(shinyjs)
library(shinyWidgets)
library(readxl)
library(corrplot)
library(igraph)
library(ggplot2)
library(dplyr)
library(viridis)

# Enable WGCNA options
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

server <- function(input, output, session) {
  
  # Enhanced reactive values to store all data and results
  values <- reactiveValues(
    expr_data = NULL,
    filtered_data = NULL,
    metadata = NULL,
    wgcna_results = NULL,
    power_results = NULL,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    current_progress = "",
    module_colors = NULL,
    gene_names = NULL,
    sample_names = NULL,
    outlier_samples = NULL,
    filtering_applied = FALSE
  )
  
  # Enhanced file upload handling with multiple formats
  observe({
    req(input$file)
    
    values$analysis_running <- TRUE
    values$current_progress <- "Loading data..."
    
    tryCatch({
      # Determine file type and load accordingly
      file_ext <- tools::file_ext(input$file$datapath)
      
      if (file_ext %in% c("csv", "txt", "tsv")) {
        # Determine separator
        sep <- switch(input$separator %||% ",",
                     "," = ",",
                     "\t" = "\t", 
                     ";" = ";")
        
        values$expr_data <- read.table(
          input$file$datapath,
          header = input$hasHeader %||% TRUE,
          row.names = if(input$hasRownames %||% TRUE) 1 else NULL,
          sep = sep,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else if (file_ext == "xlsx") {
        values$expr_data <- read_excel(input$file$datapath, sheet = 1)
        if (input$hasRownames %||% TRUE) {
          rownames(values$expr_data) <- values$expr_data[[1]]
          values$expr_data <- values$expr_data[,-1]
        }
      }
      
      # Convert to numeric matrix
      numeric_cols <- sapply(values$expr_data, is.numeric)
      if (!all(numeric_cols)) {
        values$expr_data <- values$expr_data[, numeric_cols]
        showNotification("Non-numeric columns removed", type = "warning")
      }
      
      # Store gene and sample names
      values$gene_names <- rownames(values$expr_data)
      values$sample_names <- colnames(values$expr_data)
      values$filtered_data <- values$expr_data
      
      values$analysis_running <- FALSE
      showNotification("🎉 Data loaded successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      values$analysis_running <- FALSE
      showNotification(paste("❌ Error loading file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Handle metadata file upload (currently disabled in UI)
  # observe({
  #   req(input$metaFile)
  #   tryCatch({
  #     values$metadata <- read.csv(input$metaFile$datapath, stringsAsFactors = FALSE)
  #     showNotification("📋 Metadata loaded successfully!", type = "success")
  #   }, error = function(e) {
  #     showNotification(paste("Error loading metadata:", e$message), type = "error")
  #   })
  # })
  
  # Reactive outputs for UI status
  output$analysisRunning <- reactive({
    values$analysis_running
  })
  outputOptions(output, "analysisRunning", suspendWhenHidden = FALSE)
  
  output$fileUploaded <- reactive({
    !is.null(values$expr_data)
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  output$networkAvailable <- reactive({
    values$analysis_complete && !is.null(values$wgcna_results)
  })
  outputOptions(output, "networkAvailable", suspendWhenHidden = FALSE)
  
  output$modulesAvailable <- reactive({
    values$analysis_complete && !is.null(values$module_colors)
  })
  outputOptions(output, "modulesAvailable", suspendWhenHidden = FALSE)
  
  output$resultsAvailable <- reactive({
    values$analysis_complete && !is.null(values$wgcna_results)
  })
  outputOptions(output, "resultsAvailable", suspendWhenHidden = FALSE)
  
  # Progress text output
  output$progressText <- renderText({
    values$current_progress
  })
  
  # Upload summary
  output$uploadSummary <- renderText({
    if (!is.null(values$expr_data)) {
      paste0("📊 ", nrow(values$expr_data), " genes × ", 
             ncol(values$expr_data), " samples loaded")
    }
  })
  
  # Enhanced data preview with better formatting
  output$dataPreview <- DT::renderDataTable({
    req(values$filtered_data)
    
    preview_data <- head(values$filtered_data, 100)
    preview_data <- round(preview_data, 3)
    
    DT::datatable(
      preview_data,
      options = list(
        scrollX = TRUE,
        scrollY = "350px",
        pageLength = 15,
        autoWidth = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = "_all"),
          list(width = '80px', targets = "_all")
        ),
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      class = "cell-border stripe hover",
      rownames = TRUE
    ) %>% 
    DT::formatRound(columns = 1:ncol(preview_data), digits = 3) %>%
    DT::formatStyle(
      columns = 1:ncol(preview_data),
      backgroundColor = styleInterval(
        cuts = c(0, 5, 10), 
        values = c('#ffe6e6', '#fff2e6', '#e6ffe6', '#e6f3ff')
      )
    )
  })
  
  # Additional data visualizations for upload tab
  output$expressionDist <- renderPlotly({
    req(values$filtered_data)
    
    # Sample a subset for visualization
    sample_data <- values$filtered_data[sample(nrow(values$filtered_data), 
                                              min(1000, nrow(values$filtered_data))), ]
    
    p <- ggplot(data.frame(expression = as.vector(as.matrix(sample_data))), 
                aes(x = expression)) +
      geom_histogram(bins = 50, fill = "#3c8dbc", alpha = 0.7) +
      labs(title = "Expression Distribution", x = "Expression Level", y = "Frequency") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$sampleCorrelation <- renderPlotly({
    req(values$filtered_data)
    
    if (ncol(values$filtered_data) > 50) {
      # Sample subset for large datasets
      sample_subset <- sample(ncol(values$filtered_data), 50)
      cor_data <- cor(values$filtered_data[, sample_subset], use = "complete.obs")
    } else {
      cor_data <- cor(values$filtered_data, use = "complete.obs")
    }
    
    plot_ly(z = cor_data, type = "heatmap", colorscale = "RdBu",
            hovertemplate = "Sample 1: %{x}<br>Sample 2: %{y}<br>Correlation: %{z:.3f}<extra></extra>") %>%
      layout(title = "Sample Correlation Heatmap")
  })
  
  output$missingDataPlot <- renderPlot({
    req(values$expr_data)
    
    # Calculate missing data percentage
    missing_data <- is.na(values$expr_data) | values$expr_data == 0
    missing_percent <- rowSums(missing_data) / ncol(values$expr_data) * 100
    
    ggplot(data.frame(missing_percent = missing_percent), aes(x = missing_percent)) +
      geom_histogram(bins = 30, fill = "#e74c3c", alpha = 0.7) +
      labs(title = "Missing/Zero Expression Data", 
           x = "Percentage Missing (%)", y = "Number of Genes") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Enhanced data summary value boxes
  output$nGenes <- renderValueBox({
    valueBox(
      value = if(!is.null(values$filtered_data)) {
        formatC(nrow(values$filtered_data), format = "d", big.mark = ",")
      } else "—",
      subtitle = "Genes (Filtered)",
      icon = icon("dna"),
      color = "blue"
    )
  })
  
  output$nSamples <- renderValueBox({
    valueBox(
      value = if(!is.null(values$filtered_data)) {
        formatC(ncol(values$filtered_data), format = "d", big.mark = ",")
      } else "—",
      subtitle = "Samples", 
      icon = icon("vials"),
      color = "green"
    )
  })
  
  output$dataSize <- renderValueBox({
    valueBox(
      value = if(!is.null(values$expr_data)) {
        size_mb <- format(object.size(values$expr_data), units = "MB")
        gsub(" Mb", "MB", size_mb)
      } else "—",
      subtitle = "Data Size",
      icon = icon("hdd"),
      color = "yellow"
    )
  })
  
  output$dataCompleteness <- renderValueBox({
    valueBox(
      value = if(!is.null(values$expr_data)) {
        completeness <- round((1 - sum(is.na(values$expr_data)) / 
                              (nrow(values$expr_data) * ncol(values$expr_data))) * 100, 1)
        paste0(completeness, "%")
      } else "—",
      subtitle = "Data Completeness",
      icon = icon("check-circle"),
      color = "purple"
    )
  })
  
  # Note: Removed deprecated renderDropdownMenu to fix rendering issues
  
  # Data filtering reactive
  observeEvent(input$applyFilters, {
    req(values$expr_data)
    
    # Apply filters here
    showNotification("Filters applied!", type = "message")
  })
  
  # Enhanced WGCNA analysis with detailed progress tracking
  wgcnaResults <- eventReactive(input$runAnalysis, {
    
    # Prevent multiple simultaneous executions
    if (values$analysis_running) {
      showNotification("Analysis already in progress...", type = "warning")
      return(NULL)
    }
    
    # Check if data is loaded
    if (is.null(values$expr_data)) {
      showNotification("Please upload expression data first!", type = "error")
      return(NULL)
    }
    
    values$analysis_running <- TRUE
    values$analysis_complete <- FALSE
    
    # Comprehensive progress tracking
    withProgress(message = "🚀 Executing WGCNA Pipeline...", value = 0, {
      
      tryCatch({
        setProgress(0.1, detail = "Initializing analysis...")
        values$current_progress <- "Initializing WGCNA analysis..."
        
        # Data preparation - use expr_data if filtered_data doesn't exist
        setProgress(0.15, detail = "Preparing expression data...")
        values$current_progress <- "Preparing expression data..."
        
        expr_matrix <- if (!is.null(values$filtered_data)) {
          as.matrix(values$filtered_data)
        } else {
          as.matrix(values$expr_data)
        }
        
        # Basic filtering - remove genes with too many missing values
        setProgress(0.2, detail = "Filtering genes...")
        values$current_progress <- "Filtering genes with missing data..."
        
        missing_prop <- rowMeans(is.na(expr_matrix))
        keep_genes <- missing_prop < 0.2  # Keep genes with <20% missing
        expr_matrix <- expr_matrix[keep_genes, ]
        
        # Select most variable genes (use a reasonable default)
        setProgress(0.25, detail = "Selecting variable genes...")
        values$current_progress <- "Selecting most variable genes..."
        
        gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
        n_genes <- min(5000, nrow(expr_matrix))  # Use up to 5000 most variable genes
        top_var_genes <- head(order(gene_vars, decreasing = TRUE), n_genes)
        expr_matrix <- expr_matrix[top_var_genes, ]
        
        # Transpose for WGCNA (samples as rows)
        expr_matrix <- t(expr_matrix)
      
        # Check data quality
        setProgress(0.3, detail = "Checking data quality...")
        values$current_progress <- "Performing data quality checks..."
        
        gsg <- goodSamplesGenes(expr_matrix, verbose = 0)
        if (!gsg$allOK) {
          expr_matrix <- expr_matrix[gsg$goodSamples, gsg$goodGenes]
          showNotification("Removed problematic genes/samples", type = "warning")
        }
        
        # Use the soft power from UI input
        setProgress(0.4, detail = "Using specified soft-thresholding power...")
        values$current_progress <- "Using selected soft power..."
        
        power_selected <- input$softPower
        
        setProgress(0.5, detail = paste("Using power:", power_selected))
        values$current_progress <- paste("Selected soft power:", power_selected)
        
        # Network construction
        setProgress(0.6, detail = "Constructing co-expression network...")
        values$current_progress <- "Building adjacency matrix..."
        
        # Simplified network construction
        bwnet <- blockwiseModules(
          expr_matrix,
          power = power_selected,
          networkType = input$networkType,
          deepSplit = 2,
          minModuleSize = input$minModuleSize,
          mergeCutHeight = input$mergeCutHeight,
          numericLabels = TRUE,
          saveTOMs = FALSE,
          verbose = 0
        )
        
        module_colors <- labels2colors(bwnet$colors)
        
        setProgress(0.8, detail = "Finalizing results...")
        values$current_progress <- "Organizing results..."
        
        # Store results
        results <- list(
          expr_data = expr_matrix,
          module_colors = module_colors,
          module_labels = bwnet$colors,
          dendrograms = bwnet$dendrograms,
          power_used = power_selected,
          n_modules = length(unique(bwnet$colors)) - 1  # -1 for grey module
        )
        
        values$wgcna_results <- results
        values$module_colors <- module_colors
        values$analysis_running <- FALSE
        values$analysis_complete <- TRUE
        
        setProgress(1.0, detail = "Analysis complete!")
        values$current_progress <- "WGCNA analysis completed!"
        
        showNotification("WGCNA analysis completed successfully!", 
                        type = "message", duration = 5)
        
        return(results)
        
      }, error = function(e) {
        values$analysis_running <- FALSE
        values$current_progress <- paste("Error:", e$message)
        showNotification(paste("Analysis failed:", e$message), 
                        type = "error", duration = 10)
        return(NULL)
      })
    })
  })
  
  # Quick analysis with default parameters
  quickAnalysisResults <- eventReactive(input$quickAnalysis, {
    req(values$filtered_data)
    
    values$analysis_running <- TRUE
    values$analysis_complete <- FALSE
    
    withProgress(message = "⚡ Quick WGCNA Analysis...", value = 0, {
      
      setProgress(0.2, detail = "Using default parameters...")
      values$current_progress <- "Running quick analysis with optimal defaults..."
      
      expr_matrix <- t(as.matrix(values$filtered_data))
      
      # Quick power selection
      setProgress(0.4, detail = "Auto-selecting power...")
      sft <- pickSoftThreshold(expr_matrix, verbose = 0)
      power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
      
      setProgress(0.6, detail = "Building network...")
      
      # Quick blockwise analysis
      bwnet <- blockwiseModules(
        expr_matrix,
        power = power,
        networkType = "unsigned", 
        deepSplit = 2,
        minModuleSize = 30,
        mergeCutHeight = 0.25,
        numericLabels = TRUE,
        saveTOMs = FALSE,
        verbose = 0
      )
      
      setProgress(1, detail = "Quick analysis complete!")
      values$current_progress <- "Quick analysis completed!"
      
      list(
        expr_data = expr_matrix,
        power_selected = power,
        module_colors = labels2colors(bwnet$colors),
        module_labels = bwnet$colors, 
        MEs = bwnet$MEs,
        gene_tree = bwnet$dendrograms[[1]],
        power_results = sft,
        network_type = "unsigned",
        correlation_type = "pearson"
      )
    })
  })
  
  # Observe analysis completion
  # observe({
  #   if (!is.null(wgcnaResults())) {
  #     values$wgcna_results <- wgcnaResults()
  #     values$analysis_complete <- TRUE
  #     values$analysis_running <- FALSE
  #     values$module_colors <- values$wgcna_results$module_colors
  #     
  #     showNotification("🎉 WGCNA analysis completed successfully!", 
  #                     type = "message", duration = 5)
  #   }
  # })
  
  observe({
    if (!is.null(quickAnalysisResults())) {
      values$wgcna_results <- quickAnalysisResults()
      values$analysis_complete <- TRUE
      values$analysis_running <- FALSE
      values$module_colors <- values$wgcna_results$module_colors
      
      showNotification("⚡ Quick analysis completed!", 
                      type = "success", duration = 3)
    }
  })
  
  # Full WGCNA Analysis (observeEvent with progress)
  observeEvent(input$runFullAnalysis, {
    req(values$expr_data)
    
    values$analysis_running <- TRUE
    values$analysis_complete <- FALSE
    
    # Create progress indicator
    progress <- Progress$new(session, min = 0, max = 1)
    on.exit(progress$close())
    
    tryCatch({
      # Step 1: Filter genes
      progress$set(message = "Filtering genes...", value = 0.2)
      expr_filtered <- filterGenesByVariance(values$expr_data, topN = input$topVar)
      
      # Step 2: Normalize
      progress$set(message = "Normalizing data...", value = 0.4)
      expr_norm <- normalizeExpression(expr_filtered, method = input$normMethod)
      
      # Step 3: Power selection
      progress$set(message = "Selecting soft threshold power...", value = 0.6)
      datExpr <- t(expr_norm)
      power_range <- input$powerRange[1]:input$powerRange[2]
      power_res <- pickSoftPower(
        datExpr, 
        powers = power_range,
        networkType = input$networkType,
        corType = input$corType
      )
      values$power_results <- power_res
      
      # Step 4: Network construction and module detection
      progress$set(message = "Constructing network and detecting modules...", value = 0.8)
      result <- generateRobustDendrogram(
        exprData = expr_norm,
        power = power_res$power,
        deepSplit = input$deepSplit,
        minModuleSize = input$minModSize,
        mergeCutHeight = input$mergeHeight,
        bootstrap = FALSE,
        tomType = input$networkType,
        corType = input$corType
      )
      
      progress$set(message = "Analysis complete!", value = 1.0)
      
      # Store additional info
      result$datExpr <- datExpr
      result$expr_norm <- expr_norm
      result$n_genes_original <- nrow(values$expr_data)
      result$n_genes_filtered <- nrow(expr_filtered)
      result$power_diagnostics <- power_res$diagnostics
      
      values$analysis_complete <- TRUE
      values$wgcna_results <- result
      
      showNotification("WGCNA analysis completed successfully!", type = "message", duration = 5)
      
      result
      
    }, error = function(e) {
      progress$set(message = "Analysis failed!", value = 1.0)
      showNotification(paste("Analysis failed:", e$message), type = "error", duration = 10)
      NULL
    })
  })
  
  # Check if analysis is complete
  output$analysisComplete <- reactive({
    values$analysis_complete
  })
  outputOptions(output, "analysisComplete", suspendWhenHidden = FALSE)
  
  # Power selection plot
  output$powerPlot <- renderPlotly({
    req(values$power_results)
    
    fit_data <- values$power_results$fitIndices
    powers <- fit_data[, 1]
    R2 <- fit_data[, 2] 
    mean_k <- fit_data[, 5]
    
    # Create subplot
    p1 <- plot_ly(x = powers, y = R2, type = 'scatter', mode = 'lines+markers',
                  name = 'Scale Free R²', line = list(color = '#1f77b4')) %>%
      add_trace(x = powers, y = rep(0.8, length(powers)), 
                mode = 'lines', name = 'R² = 0.8', 
                line = list(color = 'red', dash = 'dash')) %>%
      layout(title = "Scale-Free Topology Fit",
             xaxis = list(title = "Soft Threshold (power)"),
             yaxis = list(title = "Scale Free R²"))
    
    p2 <- plot_ly(x = powers, y = mean_k, type = 'scatter', mode = 'lines+markers',
                  name = 'Mean Connectivity', line = list(color = '#ff7f0e')) %>%
      layout(title = "Mean Connectivity",
             xaxis = list(title = "Soft Threshold (power)"),
             yaxis = list(title = "Mean Connectivity"))
    
    subplot(p1, p2, nrows = 2, shareX = TRUE) %>%
      layout(title = "Power Selection Diagnostics")
  })
  
  # Main dendrogram plot
  output$dendPlot <- renderPlot({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    res <- values$wgcna_results
    
    # Enhanced dendrogram plotting
    par(mar = c(5, 4, 4, 2) + 0.1)
    plotDendroAndColors(
      res$geneTree,
      res$moduleColors,
      groupLabels = "Modules",
      main = "Gene Clustering Dendrogram with Module Colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      cex.colorLabels = 0.8
    )
  }, height = 500)
  
  # Interactive dendrogram
  output$dendPlotly <- renderPlotly({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    res <- values$wgcna_results
    
    tryCatch({
      plotInteractiveDendrogram(res$geneTree, res$stability)
    }, error = function(e) {
      # Fallback plot if interactive fails
      plot_ly() %>%
        add_text(x = 0.5, y = 0.5, text = "Interactive dendrogram not available", 
                 textfont = list(size = 16, color = "gray")) %>%
        layout(title = "Interactive Dendrogram", 
               xaxis = list(visible = FALSE), 
               yaxis = list(visible = FALSE))
    })
  })
  
  # Module eigengenes plot
  output$eigengenePlot <- renderPlotly({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    res <- values$wgcna_results
    
    # Calculate module eigengenes
    MEList <- moduleEigengenes(res$datExpr, colors = res$moduleColors)
    MEs <- MEList$eigengenes
    
    if(ncol(MEs) > 1) {
      # Correlation heatmap
      ME_cor <- cor(MEs, use = "pairwise.complete.obs")
      
      plot_ly(
        z = ME_cor,
        type = "heatmap",
        colorscale = "RdYlBu",
        reversescale = TRUE,
        showscale = TRUE
      ) %>%
        layout(
          title = "Module Eigengene Correlation",
          xaxis = list(title = "Modules"),
          yaxis = list(title = "Modules")
        )
    } else {
      plot_ly() %>%
        add_text(x = 0.5, y = 0.5, text = "Need multiple modules for eigengene analysis", 
                 textfont = list(size = 16, color = "gray")) %>%
        layout(title = "Module Eigengenes", 
               xaxis = list(visible = FALSE), 
               yaxis = list(visible = FALSE))
    }
  })
  
  # Module network plot
  output$moduleNetworkPlot <- renderPlotly({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    res <- values$wgcna_results
    
    # Simple module size visualization
    module_sizes <- table(res$moduleColors)
    module_sizes <- module_sizes[names(module_sizes) != "grey"]
    
    if(length(module_sizes) > 0) {
      plot_ly(
        x = names(module_sizes),
        y = as.numeric(module_sizes),
        type = "bar",
        marker = list(color = names(module_sizes))
      ) %>%
        layout(
          title = "Module Sizes",
          xaxis = list(title = "Module"),
          yaxis = list(title = "Number of Genes")
        )
    } else {
      plot_ly() %>%
        add_text(x = 0.5, y = 0.5, text = "No modules detected", 
                 textfont = list(size = 16, color = "gray")) %>%
        layout(title = "Module Network", 
               xaxis = list(visible = FALSE), 
               yaxis = list(visible = FALSE))
    }
  })
  
  # Module summary table
  output$moduleSummary <- DT::renderDataTable({
    res <- wgcnaResults()
    req(res)
    
    module_sizes <- table(res$moduleColors)
    module_df <- data.frame(
      Module = names(module_sizes),
      Size = as.numeric(module_sizes),
      Percentage = round(100 * as.numeric(module_sizes) / length(res$moduleColors), 2),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      module_df,
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      class = "cell-border stripe"
    ) %>%
      DT::formatStyle(
        "Module",
        backgroundColor = DT::styleEqual(module_df$Module, module_df$Module)
      )
  })
  
  # Quality control metrics
  output$selectedPower <- renderText({
    res <- wgcnaResults()
    if(!is.null(res)) as.character(res$power) else "—"
  })
  
  output$scaleFreeFit <- renderText({
    if(!is.null(values$power_results)) {
      power_idx <- which(values$power_results$fitIndices[,1] == values$power_results$power)
      if(length(power_idx) > 0) {
        sprintf("%.3f", values$power_results$fitIndices[power_idx, 2])
      } else "—"
    } else "—"
  })
  
  output$numModules <- renderText({
    res <- wgcnaResults()
    if(!is.null(res)) {
      unique_modules <- unique(res$moduleColors)
      as.character(length(unique_modules) - ("grey" %in% unique_modules))
    } else "—"
  })
  
  # QC Plots
  output$sampleQCPlot <- renderPlot({
    res <- wgcnaResults()
    req(res)
    
    # Sample dendrogram for QC
    sampleTree <- hclust(dist(res$datExpr), method = "average")
    par(mar = c(4, 4, 4, 2))
    plot(sampleTree, main = "Sample Clustering", 
         xlab = "", sub = "", cex = 0.7)
  })
  
  output$geneQCPlot <- renderPlot({
    res <- wgcnaResults()
    req(res)
    
    # Gene expression distribution
    expr_means <- colMeans(res$expr_norm, na.rm = TRUE)
    hist(expr_means, breaks = 50, 
         main = "Gene Expression Distribution",
         xlab = "Mean Expression Level",
         col = "lightblue", border = "white")
  })
  
  # Analysis summary
  output$analysisSummary <- renderText({
    res <- wgcnaResults()
    if(is.null(res)) return("No analysis results available.")
    
    paste(
      "WGCNA Analysis Summary",
      "======================",
      "",
      sprintf("Original genes: %d", res$n_genes_original %||% "N/A"),
      sprintf("Filtered genes: %d", res$n_genes_filtered %||% "N/A"),
      sprintf("Samples: %d", nrow(res$datExpr) %||% "N/A"),
      sprintf("Selected power: %d", res$power %||% "N/A"),
      sprintf("Modules detected: %d", length(unique(res$moduleColors)) - ("grey" %in% res$moduleColors)),
      sprintf("Unassigned genes: %.1f%%", 100 * sum(res$moduleColors == "grey") / length(res$moduleColors)),
      "",
      "Analysis completed successfully!",
      sep = "\n"
    )
  })
  
  # Download handlers
  output$downloadModules <- downloadHandler(
    filename = function() {
      paste("wgcna_modules_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      res <- wgcnaResults()
      req(res)
      
      module_df <- data.frame(
        Gene = names(res$moduleColors),
        Module = res$moduleColors,
        stringsAsFactors = FALSE
      )
      write.csv(module_df, file, row.names = FALSE)
    }
  )
  
  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste("wgcna_network_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      res <- wgcnaResults()
      req(res)
      saveRDS(res, file)
    }
  )
  
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste("wgcna_plots_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      res <- wgcnaResults()
      req(res)
      
      pdf(file, width = 12, height = 8)
      
      # Plot 1: Dendrogram
      plotDendroAndColors(
        res$geneTree,
        res$moduleColors,
        groupLabels = "Modules",
        main = "Gene Clustering Dendrogram with Module Colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05
      )
      
      # Plot 2: Power selection
      if(!is.null(values$power_results)) {
        fit_data <- values$power_results$fitIndices
        par(mfrow = c(1, 2))
        plot(fit_data[,1], fit_data[,2], type = "b",
             xlab = "Soft Threshold (power)", ylab = "Scale Free R²",
             main = "Scale-Free Topology Fit")
        abline(h = 0.8, col = "red", lty = 2)
        
        plot(fit_data[,1], fit_data[,5], type = "b",
             xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
             main = "Mean Connectivity")
      }
      
      dev.off()
    }
  )
  
  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Data quality outputs
  output$expressionRange <- renderText({
    if (!is.null(values$filtered_data)) {
      range_vals <- range(values$filtered_data, na.rm = TRUE)
      paste0("Min: ", round(range_vals[1], 2), " | Max: ", round(range_vals[2], 2))
    }
  })
  
  output$sampleQuality <- renderText({
    if (!is.null(values$filtered_data)) {
      # Basic quality assessment
      cor_matrix <- cor(values$filtered_data, use = "complete.obs")
      mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
      paste0("Mean correlation: ", round(mean_cor, 3))
    }
  })
  
  output$normalizationStatus <- renderText({
    if (!is.null(values$filtered_data)) {
      # Check if data appears normalized
      max_val <- max(values$filtered_data, na.rm = TRUE)
      if (max_val > 50) {
        "Raw counts detected"
      } else if (max_val > 20) {
        "Log-transformed"
      } else {
        "Normalized"
      }
    }
  })
  
  # Filtering preview outputs
  output$originalGenes <- renderText({
    if (!is.null(values$expr_data)) {
      formatC(nrow(values$expr_data), format = "d", big.mark = ",")
    }
  })
  
  output$originalSamples <- renderText({
    if (!is.null(values$expr_data)) {
      formatC(ncol(values$expr_data), format = "d", big.mark = ",")
    }
  })
  
  output$filteredGenes <- renderText({
    if (!is.null(values$filtered_data)) {
      formatC(nrow(values$filtered_data), format = "d", big.mark = ",")
    }
  })
  
  output$filteredSamples <- renderText({
    if (!is.null(values$filtered_data)) {
      formatC(ncol(values$filtered_data), format = "d", big.mark = ",")
    }
  })
  
  # Real-time filtering preview
  output$genesAfterFilter <- renderText({
    req(values$expr_data)
    # Simulate filtering based on current parameters
    min_expr <- input$minExpression %||% 1
    min_samples <- input$minSamples %||% 5
    
    # Count genes that would pass filter
    expressing_samples <- rowSums(values$expr_data >= min_expr, na.rm = TRUE)
    genes_pass <- sum(expressing_samples >= min_samples)
    
    formatC(genes_pass, format = "d", big.mark = ",")
  })
  
  output$samplesAfterQC <- renderText({
    req(values$expr_data)
    # Simulate sample QC
    if (input$removeOutliers %||% TRUE) {
      # Basic outlier detection simulation
      sample_means <- colMeans(values$expr_data, na.rm = TRUE)
      outliers <- abs(scale(sample_means)) > (input$outlierThreshold %||% 3)
      samples_keep <- ncol(values$expr_data) - sum(outliers)
    } else {
      samples_keep <- ncol(values$expr_data)
    }
    
    formatC(samples_keep, format = "d", big.mark = ",")
  })
  
  # Apply filtering
  observeEvent(input$applyFilters, {
    req(values$expr_data)
    
    withProgress(message = "Applying filters...", value = 0, {
      
      setProgress(0.2, detail = "Filtering genes by expression...")
      
      # Gene filtering
      min_expr <- input$minExpressionLevel %||% input$minExpression %||% 1
      min_samples <- input$minExpressingSamples %||% input$minSamples %||% 5
      
      # Filter genes by expression
      expressing_samples <- rowSums(values$expr_data >= min_expr, na.rm = TRUE)
      genes_keep <- expressing_samples >= min_samples
      
      setProgress(0.4, detail = "Filtering by variance...")
      
      # Filter by variance
      if (input$removeLowVar %||% TRUE) {
        gene_vars <- apply(values$expr_data[genes_keep, ], 1, var, na.rm = TRUE)
        var_threshold <- quantile(gene_vars, probs = (input$varianceThreshold %||% 75) / 100, na.rm = TRUE)
        genes_keep[genes_keep] <- gene_vars >= var_threshold
      }
      
      setProgress(0.6, detail = "Sample quality control...")
      
      # Sample QC
      filtered_expr <- values$expr_data[genes_keep, ]
      samples_keep <- rep(TRUE, ncol(filtered_expr))
      
      if (input$enableSampleQC %||% input$removeOutliers %||% TRUE) {
        # Outlier detection
        sample_means <- colMeans(filtered_expr, na.rm = TRUE)
        outlier_threshold <- input$outlierSDThreshold %||% input$outlierThreshold %||% 2.5
        outliers <- abs(scale(sample_means)) > outlier_threshold
        samples_keep <- !outliers
        
        values$outlier_samples <- colnames(filtered_expr)[outliers]
      }
      
      setProgress(0.8, detail = "Finalizing filtered data...")
      
      # Apply all filters
      values$filtered_data <- filtered_expr[, samples_keep]
      values$filtering_applied <- TRUE
      
      setProgress(1, detail = "Complete!")
      
      showNotification(
        paste0("✅ Filtering complete! ", 
               sum(genes_keep), " genes and ", 
               sum(samples_keep), " samples retained."),
        type = "success", duration = 4
      )
    })
  })
  
  # Sample clustering for QC
  output$sampleClusterPlot <- renderPlotly({
    req(values$filtered_data)
    
    # Hierarchical clustering of samples
    if (ncol(values$filtered_data) > 2) {
      sample_dist <- dist(t(values$filtered_data))
      sample_hclust <- hclust(sample_dist, method = "average")
      
      dend_data <- ggdendro::dendro_data(sample_hclust)
      
      p <- ggplot(dend_data$segments) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = dend_data$labels, 
                  aes(x = x, y = y, label = label), 
                  hjust = 1, angle = 90, size = 3) +
        labs(title = "Sample Clustering Dendrogram", x = "Samples", y = "Distance") +
        theme_minimal() +
        theme(axis.text.x = element_blank())
      
      ggplotly(p)
    }
  })
  
  # Outlier detection plots
  output$outlierScatterPlot <- renderPlotly({
    req(values$filtered_data)
    
    sample_means <- colMeans(values$filtered_data, na.rm = TRUE)
    sample_sds <- apply(values$filtered_data, 2, sd, na.rm = TRUE)
    
    outlier_df <- data.frame(
      sample = colnames(values$filtered_data),
      mean_expr = sample_means,
      sd_expr = sample_sds,
      is_outlier = abs(scale(sample_means)) > 2.5
    )
    
    p <- ggplot(outlier_df, aes(x = mean_expr, y = sd_expr, color = is_outlier)) +
      geom_point(size = 2, alpha = 0.7) +
      labs(title = "Sample Quality: Mean vs SD", 
           x = "Mean Expression", y = "Standard Deviation") +
      scale_color_manual(values = c("FALSE" = "#3c8dbc", "TRUE" = "#e74c3c"),
                        name = "Outlier") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$outlierBoxPlot <- renderPlot({
    req(values$filtered_data)
    
    sample_means <- colMeans(values$filtered_data, na.rm = TRUE)
    
    ggplot(data.frame(means = sample_means), aes(y = means)) +
      geom_boxplot(fill = "#3c8dbc", alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(title = "Sample Mean Expression Distribution", y = "Mean Expression") +
      theme_minimal() +
      theme(axis.text.x = element_blank())
  })
  
  # Sample correlation heatmap
  output$sampleCorHeatmap <- renderPlotly({
    req(values$filtered_data)
    
    if (ncol(values$filtered_data) <= 100) {
      cor_matrix <- cor(values$filtered_data, use = "complete.obs")
      
      plot_ly(z = cor_matrix, type = "heatmap", 
              colorscale = list(c(0, "blue"), c(0.5, "white"), c(1, "red")),
              hovertemplate = "Sample 1: %{x}<br>Sample 2: %{y}<br>Correlation: %{z:.3f}<extra></extra>") %>%
        layout(title = "Sample-Sample Correlation Matrix")
    } else {
      # For large datasets, show a message
      plot_ly() %>%
        add_annotations(text = "Too many samples for heatmap visualization\n(>100 samples)",
                       showarrow = FALSE, font = list(size = 16))
    }
  })
  
  # ============================================================================
  # NETWORK ANALYSIS OUTPUTS
  # ============================================================================
  
  # Power selection plot
  output$powerSelectionPlot <- renderPlotly({
    req(values$power_results)
    
    sft_data <- values$power_results$fitIndices
    
    # Scale-free topology fit
    p1 <- ggplot(sft_data, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
      labs(title = "Scale-Free Topology Model Fit", 
           x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit, signed R²") +
      theme_minimal()
    
    ggplotly(p1)
  })
  
  # Network connectivity plots
  output$connectivityPlot <- renderPlotly({
    req(values$power_results)
    
    sft_data <- values$power_results$fitIndices
    
    p <- ggplot(sft_data, aes(x = Power, y = mean.k.)) +
      geom_point(size = 3, alpha = 0.7, color = "#3c8dbc") +
      geom_line(alpha = 0.5, color = "#3c8dbc") +
      labs(title = "Mean Connectivity", x = "Soft Threshold (power)", y = "Mean Connectivity") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$scaleFreePlot <- renderPlotly({
    req(values$power_results)
    
    sft_data <- values$power_results$fitIndices
    
    p <- ggplot(sft_data, aes(x = log10(mean.k.), y = log10(p.k.))) +
      geom_point(size = 3, alpha = 0.7, color = "#e74c3c") +
      geom_smooth(method = "lm", se = FALSE, color = "#2c3e50") +
      labs(title = "Scale-Free Topology Check", 
           x = "log10(k)", y = "log10(p(k))") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$meanConnectivityPlot <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$TOM)) {
      # Calculate connectivity from TOM
      connectivity <- rowSums(values$wgcna_results$TOM) - 1
      module_colors <- values$wgcna_results$module_colors
      
      conn_df <- data.frame(
        connectivity = connectivity,
        module = module_colors,
        gene = names(connectivity)
      )
      
      p <- ggplot(conn_df, aes(x = module, y = connectivity, fill = module)) +
        geom_boxplot(alpha = 0.7) +
        scale_fill_identity() +
        labs(title = "Connectivity by Module", x = "Module", y = "Connectivity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggplotly(p)
    } else {
      plot_ly() %>%
        add_annotations(text = "TOM data not available", showarrow = FALSE)
    }
  })
  
  # TOM heatmap
  output$tomHeatmap <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$TOM)) {
      # Sample subset for visualization
      n_genes <- min(500, nrow(values$wgcna_results$TOM))
      sample_indices <- sample(nrow(values$wgcna_results$TOM), n_genes)
      
      tom_subset <- values$wgcna_results$TOM[sample_indices, sample_indices]
      
      plot_ly(z = tom_subset, type = "heatmap", 
              colorscale = "Viridis",
              hovertemplate = "Gene 1: %{x}<br>Gene 2: %{y}<br>TOM: %{z:.3f}<extra></extra>") %>%
        layout(title = paste("TOM Heatmap (", n_genes, "genes sampled)"))
    } else {
      plot_ly() %>%
        add_annotations(text = "TOM heatmap not available\n(Enable for smaller datasets)", 
                       showarrow = FALSE, font = list(size = 16))
    }
  })
  
  # Network graph visualization
  output$networkGraphPlot <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$TOM)) {
      # Create network graph from top connections
      tom <- values$wgcna_results$TOM
      n_genes <- min(100, nrow(tom))
      sample_indices <- sample(nrow(tom), n_genes)
      
      tom_subset <- tom[sample_indices, sample_indices]
      threshold <- quantile(tom_subset[upper.tri(tom_subset)], 0.95)
      
      # Create adjacency matrix
      adj_matrix <- tom_subset > threshold
      
      # Create igraph object
      g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
      
      if (length(E(g)) > 0) {
        # Layout
        layout <- layout_with_fr(g)
        
        # Prepare data for plotly
        edge_shapes <- list()
        for (i in 1:length(E(g))) {
          v0 <- ends(g, E(g)[i])[1]
          v1 <- ends(g, E(g)[i])[2]
          
          edge_shapes[[i]] <- list(
            type = "line",
            line = list(color = "#999", width = 0.5),
            x0 = layout[v0, 1], y0 = layout[v0, 2],
            x1 = layout[v1, 1], y1 = layout[v1, 2]
          )
        }
        
        # Node colors by module
        node_colors <- values$wgcna_results$module_colors[sample_indices]
        
        plot_ly(x = layout[, 1], y = layout[, 2], 
                color = node_colors, colors = node_colors,
                type = "scatter", mode = "markers",
                marker = list(size = 8),
                hovertemplate = "Gene: %{text}<br>Module: %{color}<extra></extra>",
                text = rownames(tom_subset)) %>%
          layout(title = "Gene Co-expression Network",
                 showlegend = FALSE,
                 shapes = edge_shapes,
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      } else {
        plot_ly() %>%
          add_annotations(text = "No significant connections to display", 
                         showarrow = FALSE)
      }
    } else {
      plot_ly() %>%
        add_annotations(text = "Network graph requires TOM data", 
                       showarrow = FALSE)
    }
  })
  
  # ============================================================================
  # MODULE ANALYSIS OUTPUTS  
  # ============================================================================
  
  # Dendrogram plot
  output$dendrogramPlot <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$gene_tree)) {
      gene_tree <- values$wgcna_results$gene_tree
      module_colors <- values$wgcna_results$module_colors
      
      # Convert to ggdendro format
      dend_data <- ggdendro::dendro_data(gene_tree)
      
      # Create base dendrogram
      p <- ggplot(dend_data$segments) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        labs(title = "Gene Clustering Dendrogram with Module Colors",
             x = "Genes", y = "Distance") +
        theme_minimal() +
        theme(axis.text.x = element_blank())
      
      # Add module color bar
      if (length(module_colors) == length(gene_tree$order)) {
        color_df <- data.frame(
          x = 1:length(gene_tree$order),
          y = -0.1,
          color = module_colors[gene_tree$order]
        )
        
        p <- p + geom_tile(data = color_df, aes(x = x, y = y, fill = color)) +
          scale_fill_identity() +
          ylim(-0.15, max(dend_data$segments$y))
      }
      
      ggplotly(p)
    }
  })
  
  # Module color visualization
  output$moduleColorPlot <- renderPlot({
    req(values$wgcna_results)
    
    module_colors <- values$wgcna_results$module_colors
    unique_colors <- unique(module_colors)
    color_counts <- table(module_colors)
    
    # Create module size plot
    module_df <- data.frame(
      module = names(color_counts),
      size = as.numeric(color_counts)
    ) %>%
      arrange(desc(size))
    
    ggplot(module_df, aes(x = reorder(module, size), y = size, fill = module)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      scale_fill_identity() +
      labs(title = "Module Sizes", x = "Module", y = "Number of Genes") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_flip()
  })
  
  # Module statistics
  output$moduleStats <- renderText({
    req(values$wgcna_results)
    
    module_colors <- values$wgcna_results$module_colors
    n_modules <- length(unique(module_colors))
    largest_module <- max(table(module_colors))
    
    paste0("Total modules: ", n_modules, "\n",
           "Largest module: ", largest_module, " genes\n",
           "Average module size: ", round(mean(table(module_colors)), 1))
  })
  
  # Top modules table
  output$topModulesTable <- renderTable({
    req(values$wgcna_results)
    
    module_sizes <- table(values$wgcna_results$module_colors)
    top_modules <- head(sort(module_sizes, decreasing = TRUE), 5)
    
    data.frame(
      Module = names(top_modules),
      Size = as.numeric(top_modules),
      Percentage = round(as.numeric(top_modules) / length(values$wgcna_results$module_colors) * 100, 1)
    )
  })
  
  # Module eigengene plot
  output$moduleEigengenePlot <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$MEs)) {
      MEs <- values$wgcna_results$MEs
      
      # Calculate module correlations
      ME_cor <- cor(MEs, use = "complete.obs")
      
      plot_ly(z = ME_cor, type = "heatmap", 
              colorscale = list(c(0, "blue"), c(0.5, "white"), c(1, "red")),
              hovertemplate = "Module 1: %{x}<br>Module 2: %{y}<br>Correlation: %{z:.3f}<extra></extra>") %>%
        layout(title = "Module Eigengene Correlations")
    }
  })
  
  # Update module choices for gene explorer
  observe({
    if (values$analysis_complete && !is.null(values$wgcna_results)) {
      module_choices <- sort(unique(values$wgcna_results$module_colors))
      updateSelectInput(session, "selectedModule", choices = module_choices)
      updateSelectInput(session, "searchModule", choices = c("All" = "all", module_choices))
      
      # Update gene choices
      gene_choices <- rownames(values$wgcna_results$expr_data)
      updateSelectizeInput(session, "searchGenes", choices = gene_choices, server = TRUE)
    }
  })
  
  # Module network plot
  output$moduleNetworkPlot <- renderPlotly({
    req(values$wgcna_results, input$selectedModule)
    
    if (!is.null(values$wgcna_results$TOM) && input$selectedModule != "") {
      # Get genes in selected module
      module_genes <- which(values$wgcna_results$module_colors == input$selectedModule)
      
      if (length(module_genes) > 0) {
        n_genes <- min(input$topGenes %||% 20, length(module_genes))
        
        # Get top connected genes in module
        if (length(module_genes) > n_genes) {
          tom_subset <- values$wgcna_results$TOM[module_genes, module_genes]
          connectivity <- rowSums(tom_subset)
          top_genes <- head(order(connectivity, decreasing = TRUE), n_genes)
          module_genes <- module_genes[top_genes]
        }
        
        tom_module <- values$wgcna_results$TOM[module_genes, module_genes]
        threshold <- quantile(tom_module[upper.tri(tom_module)], 0.8)
        
        # Create network
        adj_matrix <- tom_module > threshold
        g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
        
        if (length(E(g)) > 0) {
          layout <- layout_with_fr(g)
          
          # Edge shapes
          edge_shapes <- list()
          for (i in 1:length(E(g))) {
            v0 <- ends(g, E(g)[i])[1]
            v1 <- ends(g, E(g)[i])[2]
            
            edge_shapes[[i]] <- list(
              type = "line",
              line = list(color = input$selectedModule, width = 1),
              x0 = layout[v0, 1], y0 = layout[v0, 2],
              x1 = layout[v1, 1], y1 = layout[v1, 2]
            )
          }
          
          plot_ly(x = layout[, 1], y = layout[, 2], 
                  type = "scatter", mode = "markers",
                  marker = list(size = 10, color = input$selectedModule),
                  hovertemplate = "Gene: %{text}<extra></extra>",
                  text = rownames(tom_module)) %>%
            layout(title = paste("Module", input$selectedModule, "Network"),
                   shapes = edge_shapes,
                   showlegend = FALSE,
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        }
      }
    }
  })
  
  # Module summary table
  output$moduleSummaryTable <- DT::renderDataTable({
    req(values$wgcna_results)
    
    module_colors <- values$wgcna_results$module_colors
    
    # Calculate module statistics
    module_summary <- data.frame(
      Module = names(table(module_colors)),
      Size = as.numeric(table(module_colors)),
      stringsAsFactors = FALSE
    )
    
    # Add connectivity info if available
    if (!is.null(values$wgcna_results$TOM)) {
      connectivity_by_module <- tapply(rowSums(values$wgcna_results$TOM), 
                                      module_colors, mean)
      module_summary$Mean_Connectivity <- round(connectivity_by_module[module_summary$Module], 2)
    }
    
    # Add eigengene info if available
    if (!is.null(values$wgcna_results$MEs)) {
      # Module preservation/quality could be calculated here
      module_summary$Quality <- "Good"  # Placeholder
    }
    
    DT::datatable(
      module_summary,
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      class = "cell-border stripe hover"
    ) %>%
    DT::formatStyle(
      "Module",
      backgroundColor = styleEqual(module_summary$Module, module_summary$Module)
    )
  })
  
  # ============================================================================
  # GENE EXPLORER OUTPUTS
  # ============================================================================
  
  # Gene information table
  output$geneInfoTable <- DT::renderDataTable({
    req(values$wgcna_results, input$searchGenes)
    
    if (length(input$searchGenes) > 0) {
      gene_indices <- match(input$searchGenes, rownames(values$wgcna_results$expr_data))
      gene_indices <- gene_indices[!is.na(gene_indices)]
      
      if (length(gene_indices) > 0) {
        gene_info <- data.frame(
          Gene = rownames(values$wgcna_results$expr_data)[gene_indices],
          Module = values$wgcna_results$module_colors[gene_indices],
          stringsAsFactors = FALSE
        )
        
        # Add connectivity if available
        if (!is.null(values$wgcna_results$TOM)) {
          connectivity <- rowSums(values$wgcna_results$TOM[gene_indices, , drop = FALSE])
          gene_info$Connectivity <- round(connectivity, 2)
        }
        
        # Add expression statistics
        expr_data <- values$wgcna_results$expr_data[gene_indices, , drop = FALSE]
        gene_info$Mean_Expression <- round(rowMeans(expr_data), 2)
        gene_info$SD_Expression <- round(apply(expr_data, 1, sd), 2)
        
        DT::datatable(
          gene_info,
          options = list(pageLength = 10, autoWidth = TRUE),
          class = "cell-border stripe hover"
        )
      }
    }
  })
  
  # Gene expression profiles
  output$geneExpressionPlot <- renderPlotly({
    req(values$wgcna_results, input$searchGenes)
    
    if (length(input$searchGenes) > 0) {
      gene_indices <- match(input$searchGenes, rownames(values$wgcna_results$expr_data))
      gene_indices <- gene_indices[!is.na(gene_indices)]
      
      if (length(gene_indices) > 0) {
        expr_data <- values$wgcna_results$expr_data[gene_indices, , drop = FALSE]
        
        # Prepare data for plotting
        plot_data <- data.frame(
          Sample = rep(colnames(expr_data), each = nrow(expr_data)),
          Gene = rep(rownames(expr_data), ncol(expr_data)),
          Expression = as.vector(t(expr_data))
        )
        
        p <- ggplot(plot_data, aes(x = Sample, y = Expression, color = Gene, group = Gene)) +
          geom_line(alpha = 0.7, size = 1) +
          geom_point(alpha = 0.5) +
          labs(title = "Gene Expression Profiles", x = "Sample", y = "Expression") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggplotly(p)
      }
    }
  })
  
  # Gene co-expression network
  output$geneNetworkPlot <- renderPlotly({
    req(values$wgcna_results, input$searchGenes)
    
    if (length(input$searchGenes) > 0 && !is.null(values$wgcna_results$TOM)) {
      gene_indices <- match(input$searchGenes, rownames(values$wgcna_results$expr_data))
      gene_indices <- gene_indices[!is.na(gene_indices)]
      
      if (length(gene_indices) > 1) {
        # Get TOM subset for selected genes
        tom_subset <- values$wgcna_results$TOM[gene_indices, gene_indices]
        
        # Create network
        threshold <- quantile(tom_subset[upper.tri(tom_subset)], 0.7)
        adj_matrix <- tom_subset > threshold
        g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
        
        if (length(E(g)) > 0) {
          layout <- layout_with_kk(g)
          
          # Edge shapes
          edge_shapes <- list()
          for (i in 1:length(E(g))) {
            v0 <- ends(g, E(g)[i])[1]
            v1 <- ends(g, E(g)[i])[2]
            
            edge_shapes[[i]] <- list(
              type = "line",
              line = list(color = "#999", width = 2),
              x0 = layout[v0, 1], y0 = layout[v0, 2],
              x1 = layout[v1, 1], y1 = layout[v1, 2]
            )
          }
          
          # Node colors by module
          node_colors <- values$wgcna_results$module_colors[gene_indices]
          
          plot_ly(x = layout[, 1], y = layout[, 2],
                  color = node_colors, colors = node_colors,
                  type = "scatter", mode = "markers",
                  marker = list(size = 15),
                  hovertemplate = "Gene: %{text}<br>Module: %{color}<extra></extra>",
                  text = rownames(tom_subset)) %>%
            layout(title = "Selected Genes Co-expression Network",
                   shapes = edge_shapes,
                   showlegend = FALSE,
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        }
      }
    }
  })
  
  # Module distribution plot
  output$moduleDistributionPlot <- renderPlot({
    req(values$wgcna_results)
    
    module_sizes <- table(values$wgcna_results$module_colors)
    
    ggplot(data.frame(size = as.numeric(module_sizes)), aes(x = size)) +
      geom_histogram(bins = 20, fill = "#3c8dbc", alpha = 0.7) +
      labs(title = "Module Size Distribution", x = "Module Size", y = "Count") +
      theme_minimal()
  })
  
  # Hub genes table
  output$hubGenesTable <- renderTable({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$TOM)) {
      connectivity <- rowSums(values$wgcna_results$TOM)
      top_hubs <- head(order(connectivity, decreasing = TRUE), 10)
      
      data.frame(
        Gene = rownames(values$wgcna_results$expr_data)[top_hubs],
        Module = values$wgcna_results$module_colors[top_hubs],
        Connectivity = round(connectivity[top_hubs], 2)
      )
    }
  })
  
  # ============================================================================
  # QUALITY CONTROL OUTPUTS
  # ============================================================================
  
  # Quality metrics
  output$selectedPower <- renderText({
    if (!is.null(values$wgcna_results)) {
      as.character(values$wgcna_results$power_selected)
    } else "—"
  })
  
  output$scaleFreeFit <- renderText({
    if (!is.null(values$power_results)) {
      power_used <- values$wgcna_results$power_selected
      sft_data <- values$power_results$fitIndices
      r_squared <- sft_data$SFT.R.sq[sft_data$Power == power_used]
      if (length(r_squared) > 0) {
        paste0(round(r_squared, 3))
      } else "—"
    } else "—"
  })
  
  output$meanConnectivity <- renderText({
    if (!is.null(values$power_results)) {
      power_used <- values$wgcna_results$power_selected
      sft_data <- values$power_results$fitIndices
      mean_k <- sft_data$mean.k.[sft_data$Power == power_used]
      if (length(mean_k) > 0) {
        round(mean_k, 1)
      } else "—"
    } else "—"
  })
  
  output$numModules <- renderText({
    if (!is.null(values$wgcna_results)) {
      length(unique(values$wgcna_results$module_colors))
    } else "—"
  })
  
  output$networkR2 <- renderText({
    if (!is.null(values$power_results)) {
      power_used <- values$wgcna_results$power_selected
      sft_data <- values$power_results$fitIndices
      r_squared <- sft_data$SFT.R.sq[sft_data$Power == power_used]
      if (length(r_squared) > 0) {
        round(r_squared, 3)
      } else "—"
    } else "—"
  })
  
  # Power summary
  output$powerSummary <- renderText({
    if (!is.null(values$wgcna_results) && !is.null(values$power_results)) {
      power_used <- values$wgcna_results$power_selected
      sft_data <- values$power_results$fitIndices
      power_row <- sft_data[sft_data$Power == power_used, ]
      
      if (nrow(power_row) > 0) {
        paste0("Selected Power: ", power_used, "\n",
               "R² = ", round(power_row$SFT.R.sq, 3), "\n",
               "Mean k = ", round(power_row$mean.k., 1), "\n",
               "Network Type: ", values$wgcna_results$network_type)
      }
    }
  })
  
  # Additional QC plots
  output$networkQualityPlot <- renderPlotly({
    req(values$power_results)
    
    sft_data <- values$power_results$fitIndices
    
    p <- ggplot(sft_data, aes(x = Power)) +
      geom_line(aes(y = SFT.R.sq, color = "R²"), size = 1) +
      geom_line(aes(y = mean.k./max(mean.k.), color = "Connectivity"), size = 1) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      labs(title = "Network Quality Metrics", x = "Soft Threshold Power", y = "Value") +
      scale_color_manual(values = c("R²" = "#3c8dbc", "Connectivity" = "#e74c3c")) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$moduleStabilityPlot <- renderPlotly({
    req(values$wgcna_results)
    
    if (!is.null(values$wgcna_results$MEs)) {
      # Module eigengene variance (stability indicator)
      me_vars <- apply(values$wgcna_results$MEs, 2, var)
      
      stability_df <- data.frame(
        Module = gsub("ME", "", names(me_vars)),
        Variance = me_vars
      )
      
      p <- ggplot(stability_df, aes(x = reorder(Module, Variance), y = Variance, fill = Module)) +
        geom_bar(stat = "identity", alpha = 0.7) +
        scale_fill_identity() +
        labs(title = "Module Stability (Eigengene Variance)", x = "Module", y = "Variance") +
        theme_minimal() +
        coord_flip()
      
      ggplotly(p)
    }
  })
  
  # Sample and gene QC plots
  output$sampleQCPlot <- renderPlot({
    req(values$filtered_data)
    
    # Sample-wise statistics
    sample_means <- colMeans(values$filtered_data, na.rm = TRUE)
    sample_sds <- apply(values$filtered_data, 2, sd, na.rm = TRUE)
    
    qc_df <- data.frame(
      Sample = colnames(values$filtered_data),
      Mean = sample_means,
      SD = sample_sds
    )
    
    ggplot(qc_df, aes(x = Mean, y = SD)) +
      geom_point(alpha = 0.6, size = 2) +
      labs(title = "Sample QC: Mean vs Standard Deviation", 
           x = "Mean Expression", y = "Standard Deviation") +
      theme_minimal()
  })
  
  output$sampleClusteringQC <- renderPlot({
    req(values$filtered_data)
    
    if (ncol(values$filtered_data) > 2) {
      # Sample clustering for QC
      sample_cor <- cor(values$filtered_data, use = "complete.obs")
      sample_dist <- as.dist(1 - sample_cor)
      sample_hclust <- hclust(sample_dist, method = "average")
      
      plot(sample_hclust, main = "Sample Clustering for Quality Control",
           xlab = "Samples", sub = "", cex = 0.8)
    }
  })
  
  output$geneQCPlot <- renderPlot({
    req(values$filtered_data)
    
    gene_means <- rowMeans(values$filtered_data, na.rm = TRUE)
    gene_vars <- apply(values$filtered_data, 1, var, na.rm = TRUE)
    
    qc_df <- data.frame(Mean = gene_means, Variance = gene_vars)
    
    ggplot(qc_df, aes(x = Mean, y = Variance)) +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      labs(title = "Gene QC: Mean vs Variance", 
           x = "Mean Expression", y = "Variance") +
      theme_minimal()
  })
  
  output$geneVariancePlot <- renderPlot({
    req(values$filtered_data)
    
    gene_vars <- apply(values$filtered_data, 1, var, na.rm = TRUE)
    
    ggplot(data.frame(variance = gene_vars), aes(x = variance)) +
      geom_histogram(bins = 50, fill = "#3c8dbc", alpha = 0.7) +
      labs(title = "Gene Variance Distribution", x = "Variance", y = "Count") +
      theme_minimal()
  })
  
  # ============================================================================
  # COMPARISON AND ANALYSIS SUMMARY
  # ============================================================================
  
  output$moduleOverlapPlot <- renderPlotly({
    # Placeholder for module overlap analysis
    plot_ly() %>%
      add_annotations(text = "Module overlap analysis\ncoming soon!", 
                     showarrow = FALSE, font = list(size = 16))
  })
  
  output$parameterSensitivityPlot <- renderPlotly({
    # Placeholder for parameter sensitivity analysis
    plot_ly() %>%
      add_annotations(text = "Parameter sensitivity analysis\ncoming soon!", 
                     showarrow = FALSE, font = list(size = 16))
  })
  
  # Analysis summary
  output$analysisSummary <- renderText({
    if (values$analysis_complete && !is.null(values$wgcna_results)) {
      paste0(
        "=== WGCNA Analysis Summary ===\n\n",
        "Data: ", nrow(values$wgcna_results$expr_data), " samples × ", 
        ncol(values$wgcna_results$expr_data), " genes\n",
        "Power: ", values$wgcna_results$power_selected, "\n",
        "Network Type: ", values$wgcna_results$network_type, "\n",
        "Correlation: ", values$wgcna_results$correlation_type, "\n",
        "Modules Detected: ", length(unique(values$wgcna_results$module_colors)), "\n",
        "Largest Module: ", max(table(values$wgcna_results$module_colors)), " genes\n",
        "Analysis Date: ", Sys.time()
      )
    }
  })
  
  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================
  
  # Module assignments download
  output$downloadModules <- downloadHandler(
    filename = function() {
      paste0("WGCNA_modules_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$wgcna_results)
      
      module_data <- data.frame(
        Gene = rownames(values$wgcna_results$expr_data),
        Module = values$wgcna_results$module_colors,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(values$wgcna_results$TOM)) {
        connectivity <- rowSums(values$wgcna_results$TOM)
        module_data$Connectivity <- connectivity
      }
      
      write.csv(module_data, file, row.names = FALSE)
    }
  )
  
  # Network data download
  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste0("WGCNA_network_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(values$wgcna_results)
      saveRDS(values$wgcna_results, file)
    }
  )
  
  # Filtered expression data download
  output$downloadExpression <- downloadHandler(
    filename = function() {
      paste0("filtered_expression_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$filtered_data)
      write.csv(values$filtered_data, file)
    }
  )
  
  # TOM matrix download
  output$downloadTOM <- downloadHandler(
    filename = function() {
      paste0("TOM_matrix_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(values$wgcna_results$TOM)
      saveRDS(values$wgcna_results$TOM, file)
    }
  )
  
  # Module eigengenes download
  output$downloadEigengenes <- downloadHandler(
    filename = function() {
      paste0("module_eigengenes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$wgcna_results$MEs)
      write.csv(values$wgcna_results$MEs, file)
    }
  )
  
  # Individual plot downloads
  output$downloadDendrogram <- downloadHandler(
    filename = function() {
      paste0("dendrogram_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(values$wgcna_results)
      
      pdf(file, width = 12, height = 8)
      if (!is.null(values$wgcna_results$gene_tree)) {
        plotDendroAndColors(values$wgcna_results$gene_tree, 
                           values$wgcna_results$module_colors,
                           "Module colors", dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05)
      }
      dev.off()
    }
  )
  
  # All plots PDF download
  output$downloadAllPlots <- downloadHandler(
    filename = function() {
      paste0("WGCNA_all_plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(values$wgcna_results)
      
      pdf(file, width = 12, height = 8)
      
      # Power selection plot
      if (!is.null(values$power_results)) {
        sft_data <- values$power_results$fitIndices
        par(mfrow = c(1, 2))
        plot(sft_data$Power, -sign(sft_data$slope) * sft_data$SFT.R.sq,
             xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
             type = "n", main = "Scale-free topology")
        text(sft_data$Power, -sign(sft_data$slope) * sft_data$SFT.R.sq,
             labels = sft_data$Power, cex = 0.9, col = "red")
        abline(h = 0.80, col = "red")
        
        plot(sft_data$Power, sft_data$mean.k.,
             xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
             main = "Mean connectivity")
        text(sft_data$Power, sft_data$mean.k., labels = sft_data$Power, cex = 0.9, col = "red")
      }
      
      # Dendrogram
      if (!is.null(values$wgcna_results$gene_tree)) {
        plotDendroAndColors(values$wgcna_results$gene_tree,
                           values$wgcna_results$module_colors,
                           "Module colors", dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05,
                           main = "Gene Dendrogram and Module Colors")
      }
      
      # Module eigengene heatmap
      if (!is.null(values$wgcna_results$MEs)) {
        MEDiss <- 1 - cor(values$wgcna_results$MEs)
        METree <- hclust(as.dist(MEDiss), method = "average")
        plot(METree, main = "Module Eigengene Clustering", xlab = "", sub = "")
      }
      
      dev.off()
    }
  )
  
  # Analysis settings download
  output$downloadSettings <- downloadHandler(
    filename = function() {
      paste0("analysis_settings_", Sys.Date(), ".txt")
    },
    content = function(file) {
      settings <- list(
        power_range = input$powerRange,
        network_type = input$networkType,
        correlation_type = input$corType,
        deep_split = input$deepSplit,
        min_module_size = input$minModSize,
        merge_height = input$mergeHeight,
        normalization = input$normMethod,
        top_var_genes = input$topVar,
        analysis_date = Sys.time()
      )
      
      writeLines(capture.output(str(settings)), file)
    }
  )
  
  # Full analysis report
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste0("WGCNA_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      # This would generate a comprehensive HTML report
      # For now, create a simple summary
      req(values$wgcna_results)
      
      html_content <- paste0(
        "<html><head><title>WGCNA Analysis Report</title></head><body>",
        "<h1>WGCNA Analysis Report</h1>",
        "<h2>Analysis Summary</h2>",
        "<p>Generated on: ", Sys.time(), "</p>",
        "<p>Data: ", nrow(values$wgcna_results$expr_data), " samples × ", 
        ncol(values$wgcna_results$expr_data), " genes</p>",
        "<p>Selected Power: ", values$wgcna_results$power_selected, "</p>",
        "<p>Network Type: ", values$wgcna_results$network_type, "</p>",
        "<p>Modules Detected: ", length(unique(values$wgcna_results$module_colors)), "</p>",
        "</body></html>"
      )
      
      writeLines(html_content, file)
    }
  )
  
  # Session info download
  output$downloadSessionInfo <- downloadHandler(
    filename = function() {
      paste0("session_info_", Sys.Date(), ".txt")
    },
    content = function(file) {
      writeLines(capture.output(sessionInfo()), file)
    }
  )
}
