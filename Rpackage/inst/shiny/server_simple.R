# Simplified Server for WGCNA Explorer
# =====================================

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
  
  # Simplified reactive values
  values <- reactiveValues(
    expr_data = NULL,
    filtered_data = NULL,
    wgcna_results = NULL,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    current_progress = "",
    module_colors = NULL
  )
  
  # File upload handling
  observe({
    req(input$file)
    
    values$analysis_running <- TRUE
    values$current_progress <- "Loading data..."
    
    tryCatch({
      file_ext <- tools::file_ext(input$file$datapath)
      
      if (file_ext %in% c("csv", "txt", "tsv")) {
        values$expr_data <- read.table(
          input$file$datapath,
          header = TRUE,
          row.names = 1,
          sep = ",",
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else if (file_ext == "xlsx") {
        values$expr_data <- read_excel(input$file$datapath, sheet = 1)
        rownames(values$expr_data) <- values$expr_data[[1]]
        values$expr_data <- values$expr_data[,-1]
      }
      
      # Convert to numeric matrix
      numeric_cols <- sapply(values$expr_data, is.numeric)
      if (!all(numeric_cols)) {
        values$expr_data <- values$expr_data[, numeric_cols]
        showNotification("Non-numeric columns removed", type = "warning")
      }
      
      values$filtered_data <- values$expr_data
      values$analysis_running <- FALSE
      showNotification("🎉 Data loaded successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      values$analysis_running <- FALSE
      showNotification(paste("❌ Error loading file:", e$message), type = "error", duration = 5)
    })
  })
  
  # Reactive outputs for UI status
  output$analysisRunning <- reactive({ values$analysis_running })
  outputOptions(output, "analysisRunning", suspendWhenHidden = FALSE)
  
  output$fileUploaded <- reactive({ !is.null(values$expr_data) })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  output$networkAvailable <- reactive({ values$analysis_complete })
  outputOptions(output, "networkAvailable", suspendWhenHidden = FALSE)
  
  output$modulesAvailable <- reactive({ values$analysis_complete })
  outputOptions(output, "modulesAvailable", suspendWhenHidden = FALSE)
  
  output$resultsAvailable <- reactive({ values$analysis_complete })
  outputOptions(output, "resultsAvailable", suspendWhenHidden = FALSE)
  
  # Progress text output
  output$progressText <- renderText({ values$current_progress })
  
  # Upload summary
  output$uploadSummary <- renderText({
    if (!is.null(values$expr_data)) {
      paste0("📊 ", nrow(values$expr_data), " genes × ", 
             ncol(values$expr_data), " samples loaded")
    }
  })
  
  # Data preview
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
        autoWidth = TRUE
      ),
      class = "cell-border stripe hover",
      rownames = TRUE
    )
  })
  
  # Value boxes
  output$nGenes <- renderValueBox({
    valueBox(
      value = if(!is.null(values$filtered_data)) nrow(values$filtered_data) else "—",
      subtitle = "Genes",
      icon = icon("dna"),
      color = "blue"
    )
  })
  
  output$nSamples <- renderValueBox({
    valueBox(
      value = if(!is.null(values$filtered_data)) ncol(values$filtered_data) else "—",
      subtitle = "Samples", 
      icon = icon("vials"),
      color = "green"
    )
  })
  
  output$nModules <- renderValueBox({
    valueBox(
      value = if(values$analysis_complete && !is.null(values$module_colors)) {
        length(unique(values$module_colors)) - 1  # -1 for grey module
      } else "—",
      subtitle = "Modules",
      icon = icon("project-diagram"),
      color = "purple"
    )
  })
  
  # WGCNA Analysis - SIMPLIFIED
  observeEvent(input$runAnalysis, {
    
    # Prevent multiple executions
    if (values$analysis_running) {
      showNotification("Analysis already in progress...", type = "warning")
      return()
    }
    
    # Check if data is loaded
    if (is.null(values$expr_data)) {
      showNotification("Please upload expression data first!", type = "error")
      return()
    }
    
    values$analysis_running <- TRUE
    values$analysis_complete <- FALSE
    
    # Run analysis with progress
    withProgress(message = "🚀 Running WGCNA Analysis...", value = 0, {
      
      tryCatch({
        setProgress(0.1, detail = "Preparing data...")
        values$current_progress <- "Preparing expression data..."
        
        # Prepare data
        expr_matrix <- as.matrix(values$filtered_data)
        
        # Basic filtering
        setProgress(0.2, detail = "Filtering genes...")
        missing_prop <- rowMeans(is.na(expr_matrix))
        keep_genes <- missing_prop < 0.2
        expr_matrix <- expr_matrix[keep_genes, ]
        
        # Select most variable genes
        setProgress(0.3, detail = "Selecting variable genes...")
        gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
        n_genes <- min(1000, nrow(expr_matrix))  # Use fewer genes for speed
        top_var_genes <- head(order(gene_vars, decreasing = TRUE), n_genes)
        expr_matrix <- expr_matrix[top_var_genes, ]
        
        # Transpose for WGCNA
        expr_matrix <- t(expr_matrix)
        
        setProgress(0.4, detail = "Checking data quality...")
        gsg <- goodSamplesGenes(expr_matrix, verbose = 0)
        if (!gsg$allOK) {
          expr_matrix <- expr_matrix[gsg$goodSamples, gsg$goodGenes]
        }
        
        setProgress(0.6, detail = "Building network...")
        values$current_progress <- "Constructing co-expression network..."
        
        # Simplified network construction using blockwiseModules
        bwnet <- blockwiseModules(
          expr_matrix,
          power = input$softPower,
          networkType = input$networkType,
          deepSplit = 2,
          minModuleSize = input$minModuleSize,
          mergeCutHeight = input$mergeCutHeight,
          numericLabels = TRUE,
          saveTOMs = FALSE,
          verbose = 0
        )
        
        setProgress(0.9, detail = "Organizing results...")
        module_colors <- labels2colors(bwnet$colors)
        
        # Store results
        values$wgcna_results <- list(
          expr_data = expr_matrix,
          module_colors = module_colors,
          module_labels = bwnet$colors,
          dendrograms = bwnet$dendrograms,
          power_used = input$softPower,
          n_modules = length(unique(bwnet$colors)) - 1
        )
        
        values$module_colors <- module_colors
        values$analysis_running <- FALSE
        values$analysis_complete <- TRUE
        
        setProgress(1.0, detail = "Complete!")
        values$current_progress <- "WGCNA analysis completed!"
        
        showNotification("WGCNA analysis completed successfully!", 
                        type = "message", duration = 5)
        
      }, error = function(e) {
        values$analysis_running <- FALSE
        values$current_progress <- paste("Error:", e$message)
        showNotification(paste("Analysis failed:", e$message), 
                        type = "error", duration = 10)
      })
    })
  })
  
  # Simple dendrogram plot - only render when analysis is complete
  output$dendPlot <- renderPlot({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    res <- values$wgcna_results
    
    if (!is.null(res$dendrograms) && length(res$dendrograms) > 0) {
      par(mar = c(5, 4, 4, 2) + 0.1)
      plotDendroAndColors(
        res$dendrograms[[1]],
        res$module_colors,
        "Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05,
        cex.colorLabels = 0.8
      )
    }
  }, height = 500)
  
  # Simple module summary
  output$moduleSummary <- renderText({
    req(values$analysis_complete)
    req(values$wgcna_results)
    
    n_modules <- values$wgcna_results$n_modules
    paste0("Analysis identified ", n_modules, " co-expression modules using soft power ", 
           values$wgcna_results$power_used)
  })
  
  # Module sizes table
  output$moduleTable <- DT::renderDataTable({
    req(values$analysis_complete)
    req(values$module_colors)
    
    module_sizes <- table(values$module_colors)
    module_df <- data.frame(
      Module = names(module_sizes),
      Size = as.numeric(module_sizes),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      module_df,
      options = list(pageLength = 15, autoWidth = TRUE),
      rownames = FALSE
    )
  })
  
}
