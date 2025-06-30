# Enhanced WGCNA Explorer UI
# ============================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

# Define Enhanced UI
ui <- dashboardPage(
  skin = "blue",
  
  # Enhanced Header
  dashboardHeader(
    title = "🧬 Enhanced WGCNA Explorer",
    titleWidth = 350,
    
    # Help dropdown menu
    dropdownMenu(
      type = "notifications",
      icon = icon("question-circle"),
      badgeStatus = "info",
      headerText = "Quick Help",
      notificationItem(
        text = "Upload expression data to start analysis",
        icon = icon("upload"),
        status = "info"
      ),
      notificationItem(
        text = "Adjust parameters for optimal results", 
        icon = icon("sliders-h"),
        status = "warning"
      ),
      notificationItem(
        text = "Export results when analysis complete",
        icon = icon("download"),
        status = "success"
      )
    ),
    
    # GitHub link
    tags$li(class = "dropdown",
            tags$a(href = "https://github.com/yourusername/MyWGCNAResearchProject", 
                   target = "_blank",
                   tags$i(class = "fa fa-github", style = "font-size: 18px;"),
                   style = "padding: 15px; color: white;")
    )
  ),
  
  # Enhanced Sidebar with advanced features
  dashboardSidebar(
    width = 350,
    
    # Analysis progress indicator
    conditionalPanel(
      condition = "output.analysisRunning",
      div(style = "padding: 15px; background: #f39c12; color: white; text-align: center; margin: 10px;",
          icon("spinner", class = "fa-spin"),
          br(),
          "Analysis in Progress...",
          br(),
          textOutput("progressText")
      )
    ),
    
    sidebarMenu(
      id = "tabs",
      menuItem("📊 Data Management", tabName = "upload", icon = icon("database"),
               menuSubItem("Upload Data", tabName = "upload"),
               menuSubItem("Data Filtering", tabName = "filtering"),
               menuSubItem("Sample QC", tabName = "sample_qc")
      ),
      menuItem("⚙️ Analysis Pipeline", tabName = "settings", icon = icon("cogs"),
               menuSubItem("Preprocessing", tabName = "preprocessing"),
               menuSubItem("Network Settings", tabName = "network_settings"),
               menuSubItem("Module Detection", tabName = "module_settings")
      ),
      menuItem("🌳 Network Explorer", tabName = "network", icon = icon("project-diagram")),
      menuItem("🎨 Module Analysis", tabName = "modules", icon = icon("palette")),
      menuItem("🔍 Gene Explorer", tabName = "genes", icon = icon("search")),
      menuItem("📈 Quality Metrics", tabName = "qc", icon = icon("chart-line")),
      menuItem("� Comparisons", tabName = "compare", icon = icon("balance-scale")),
      menuItem("📋 Export Center", tabName = "export", icon = icon("download")),
      menuItem("ℹ️ Documentation", tabName = "help", icon = icon("book"))
    ),
    
    # Advanced Data Upload Panel
    conditionalPanel(
      condition = "input.tabs == 'upload'",
      br(),
      div(style = "padding: 15px;",
          h4("📁 Expression Data", style = "color: #3c8dbc;"),
          
          fileInput("exprFile", 
                   label = "Upload Expression Matrix",
                   accept = c(".csv", ".txt", ".tsv", ".xlsx"),
                   buttonLabel = "Browse...",
                   placeholder = "Select file..."),
          
          helpText("📝 Supported: CSV, TXT, TSV, Excel"),
          helpText("🔢 Format: Genes (rows) × Samples (columns)"),
          
          # Optional metadata upload
          fileInput("metaFile",
                   label = "Sample Metadata (Optional)",
                   accept = c(".csv", ".txt"),
                   buttonLabel = "Browse...",
                   placeholder = "Optional metadata..."),
          
          # File format options
          conditionalPanel(
            condition = "input.exprFile",
            h5("📋 Import Options", style = "color: #666;"),
            checkboxInput("hasHeader", "File has header row", value = TRUE),
            checkboxInput("hasRownames", "First column is gene names", value = TRUE),
            selectInput("separator", "Field separator:",
                       choices = list("Comma (,)" = ",", "Tab" = "\t", "Semicolon (;)" = ";"),
                       selected = ",")
          ),
          
          conditionalPanel(
            condition = "output.fileUploaded",
            div(class = "alert alert-success", 
                icon("check-circle"), " Data loaded successfully!",
                br(),
                textOutput("uploadSummary")
            )
          )
      )
    ),
    
    # Data Filtering Panel
    conditionalPanel(
      condition = "input.tabs == 'filtering'",
      br(),
      div(style = "padding: 15px;",
          h4("🔍 Data Filtering", style = "color: #3c8dbc;"),
          
          conditionalPanel(
            condition = "output.fileUploaded",
            
            h5("Gene Filtering", style = "color: #666;"),
            sliderInput("minExpression", 
                       "Minimum Expression Level:",
                       min = 0, max = 20, value = 1, step = 0.1),
            
            sliderInput("minSamples", 
                       "Min Samples Expressing:",
                       min = 1, max = 50, value = 5, step = 1),
            
            sliderInput("variancePercentile", 
                       "Variance Percentile Cutoff:",
                       min = 0, max = 100, value = 75, step = 5),
            
            hr(),
            h5("Sample Filtering", style = "color: #666;"),
            checkboxInput("removeOutliers", "Auto-remove outlier samples", value = TRUE),
            
            conditionalPanel(
              condition = "input.removeOutliers",
              sliderInput("outlierThreshold", 
                         "Outlier Z-score threshold:",
                         min = 1, max = 5, value = 3, step = 0.5)
            ),
            
            actionBttn("applyFilters", 
                      "Apply Filters",
                      style = "gradient",
                      color = "primary",
                      size = "sm")
          ),
          
          conditionalPanel(
            condition = "!output.fileUploaded",
            div(class = "alert alert-info",
                icon("info-circle"), " Upload data first to enable filtering options")
          )
      )
    ),
    
    # Advanced Analysis Settings Panel
    conditionalPanel(
      condition = "input.tabs == 'preprocessing' || input.tabs == 'settings'",
      br(),
      div(style = "padding: 15px;",
          h4("🔧 Data Preprocessing", style = "color: #3c8dbc;"),
          
          pickerInput("normMethod", 
                     "Normalization Method:",
                     choices = list(
                       "Log2 Transform" = "log2",
                       "Variance Stabilizing" = "vst",
                       "Quantile Normalization" = "quantile",
                       "TMM (edgeR)" = "tmm",
                       "None" = "none"
                     ),
                     selected = "log2",
                     options = pickerOptions(
                       actionsBox = TRUE,
                       liveSearch = TRUE
                     )),
          
          # Advanced normalization options
          conditionalPanel(
            condition = "input.normMethod != 'none'",
            checkboxInput("centerScale", "Center and scale genes", value = TRUE),
            checkboxInput("removeBatch", "Batch effect correction", value = FALSE)
          ),
          
          conditionalPanel(
            condition = "input.removeBatch",
            textInput("batchVariable", "Batch column name:", placeholder = "batch")
          ),
          
          sliderInput("topVar", 
                     "Most Variable Genes:",
                     min = 500, max = 25000, value = 10000, step = 500,
                     animate = animationOptions(interval = 200)),
          
          # Real-time filtering preview
          conditionalPanel(
            condition = "output.fileUploaded",
            h5("📊 Filter Preview", style = "color: #666;"),
            div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px;",
                "Genes after filtering: ", textOutput("genesAfterFilter", inline = TRUE), br(),
                "Samples after QC: ", textOutput("samplesAfterQC", inline = TRUE)
            )
          )
      )
    ),
    
    conditionalPanel(
      condition = "input.tabs == 'network_settings'",
      br(),
      div(style = "padding: 15px;",
          h4("🌐 Network Construction", style = "color: #3c8dbc;"),
          
          # Power selection with preview
          h5("Soft-thresholding Power", style = "color: #666;"),
          sliderInput("powerRange", 
                     "Power Testing Range:",
                     min = 1, max = 35, value = c(1, 25), step = 1),
          
          pickerInput("networkType",
                     "Network Type:",
                     choices = list(
                       "Unsigned" = "unsigned",
                       "Signed" = "signed", 
                       "Signed Hybrid" = "signed hybrid"
                     ),
                     selected = "unsigned"),
          
          pickerInput("corType",
                     "Correlation Method:",
                     choices = list(
                       "Pearson" = "pearson",
                       "Biweight Midcorrelation" = "bicor",
                       "Spearman" = "spearman"
                     ),
                     selected = "pearson"),
          
          # Advanced network options
          h5("Advanced Options", style = "color: #666;"),
          checkboxInput("useBlockwise", "Use blockwise construction", value = TRUE),
          
          conditionalPanel(
            condition = "input.useBlockwise",
            numericInput("maxBlockSize", "Max block size:", value = 5000, min = 1000, max = 20000, step = 1000),
            checkboxInput("randomSeed", "Set random seed", value = TRUE),
            conditionalPanel(
              condition = "input.randomSeed",
              numericInput("seedValue", "Seed value:", value = 12345, min = 1, max = 999999)
            )
          ),
          
          # Network quality thresholds
          h5("Quality Thresholds", style = "color: #666;"),
          sliderInput("minR2", "Minimum Scale-Free R²:", min = 0.5, max = 1, value = 0.8, step = 0.05),
          sliderInput("maxMeanK", "Maximum Mean Connectivity:", min = 50, max = 500, value = 200, step = 25)
      )
    ),
    
    # Module Detection Settings Panel
    conditionalPanel(
      condition = "input.tabs == 'module_settings'",
      br(),
      div(style = "padding: 15px;",
          h4("🌳 Module Detection", style = "color: #3c8dbc;"),
          
          # Clustering method selection
          pickerInput("clusterMethod", 
                     "Clustering Method:",
                     choices = list(
                       "Hierarchical (average)" = "average",
                       "Hierarchical (complete)" = "complete", 
                       "Hierarchical (ward)" = "ward.D2"
                     ),
                     selected = "average"),
          
          sliderInput("deepSplit", 
                     "Tree Cut Sensitivity:",
                     min = 0, max = 4, value = 2, step = 1),
          
          numericInput("minModSize", 
                      "Minimum Module Size:",
                      value = 30, min = 10, max = 500, step = 5),
          
          sliderInput("mergeHeight", 
                     "Module Merge Threshold:",
                     min = 0, max = 1, value = 0.25, step = 0.02),
          
          # Advanced module options
          h5("Advanced Options", style = "color: #666;"),
          checkboxInput("detectCutHeight", "Auto-detect cut height", value = TRUE),
          checkboxInput("pamStage", "Enable PAM stage", value = TRUE),
          checkboxInput("reassignThreshold", "Reassign threshold", value = TRUE),
          
          conditionalPanel(
            condition = "input.reassignThreshold",
            sliderInput("reassignThres", "P-value threshold:", min = 0.01, max = 0.2, value = 0.05, step = 0.01)
          ),
          
          hr(),
          # Run analysis button
          div(style = "text-align: center;",
              actionBttn("runBtn", 
                        "🚀 Execute WGCNA Pipeline",
                        style = "gradient",
                        color = "primary",
                        size = "lg",
                        block = TRUE),
              br(), br(),
              
              # Quick analysis button
              actionBttn("quickAnalysis", 
                        "⚡ Quick Analysis (Default Settings)",
                        style = "gradient",
                        color = "success", 
                        size = "md",
                        block = TRUE)
          )
      )
    )
  ),
  
  # Right sidebar for advanced controls
  rightsidebar = rightSidebar(
    background = "light",
    width = 350,
    
    h4("🎛️ Advanced Controls", style = "padding: 15px; color: #3c8dbc;"),
    
    # Theme selector
    div(style = "padding: 15px;",
        h5("🎨 Theme", style = "color: #666;"),
        pickerInput("uiTheme", "UI Theme:",
                   choices = list(
                     "Blue" = "blue",
                     "Black" = "black", 
                     "Purple" = "purple",
                     "Green" = "green",
                     "Red" = "red",
                     "Yellow" = "yellow"
                   ),
                   selected = "blue"),
        
        # Plot settings
        h5("📊 Plot Settings", style = "color: #666;"),
        sliderInput("plotHeight", "Default plot height:", min = 300, max = 800, value = 500, step = 50),
        sliderInput("plotDPI", "Plot resolution (DPI):", min = 72, max = 300, value = 150, step = 25),
        
        # Export settings
        h5("💾 Export Settings", style = "color: #666;"),
        pickerInput("exportFormat", "Default format:",
                   choices = list("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"),
                   selected = "pdf"),
        
        checkboxInput("includeData", "Include raw data in exports", value = TRUE),
        checkboxInput("includeSettings", "Include analysis settings", value = TRUE)
    )
  ),

  # Enhanced Dashboard Body with Premium Features
  dashboardBody(
    # Enable shinyjs
    useShinyjs(),
    
    # Premium CSS styling
    tags$head(
      tags$style(HTML("
        .main-header .navbar {
          background-color: #367fa9 !important;
        }
        .skin-blue-light .main-header .logo {
          background-color: #367fa9 !important;
        }
        .skin-blue-light .main-header .logo:hover {
          background-color: #2e6da4 !important;
        }
        .content-wrapper, .right-side {
          background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        }
        .box {
          border-radius: 12px;
          box-shadow: 0 4px 15px rgba(0,0,0,0.1);
          transition: all 0.3s ease;
        }
        .box:hover {
          transform: translateY(-2px);
          box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }
        .nav-tabs-custom > .nav-tabs > li.active {
          border-top-color: #3c8dbc;
        }
        .progress-bar {
          background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
        }
        .alert-info {
          background: linear-gradient(45deg, #d9edf7 0%, #b8e0f5 100%);
          border-color: #bce8f1;
          color: #31708f;
          border-radius: 8px;
        }
        .metric-box {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          padding: 25px;
          border-radius: 15px;
          text-align: center;
          margin: 15px 0;
          box-shadow: 0 6px 20px rgba(102, 126, 234, 0.3);
          transition: all 0.3s ease;
        }
        .metric-box:hover {
          transform: translateY(-3px);
          box-shadow: 0 10px 30px rgba(102, 126, 234, 0.4);
        }
        .metric-value {
          font-size: 2.5em;
          font-weight: bold;
          text-shadow: 0 2px 4px rgba(0,0,0,0.3);
        }
        .metric-label {
          font-size: 1em;
          opacity: 0.95;
          margin-top: 5px;
        }
        .status-card {
          background: white;
          border-radius: 10px;
          padding: 20px;
          margin: 10px 0;
          box-shadow: 0 3px 10px rgba(0,0,0,0.1);
        }
        .analysis-progress {
          background: linear-gradient(45deg, #f39c12 0%, #e67e22 100%);
          color: white;
          padding: 20px;
          border-radius: 10px;
          text-align: center;
          margin: 15px 0;
        }
        .btn-gradient-primary {
          background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
          border: none;
          border-radius: 25px;
          padding: 12px 30px;
          color: white;
          font-weight: bold;
          transition: all 0.3s ease;
        }
        .btn-gradient-primary:hover {
          transform: translateY(-2px);
          box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
        }
        .custom-tab-content {
          background: white;
          border-radius: 12px;
          padding: 20px;
          margin: 10px 0;
          box-shadow: 0 3px 15px rgba(0,0,0,0.08);
        }
      "))
    ),

    tabItems(
      # Enhanced Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "📊 Data Overview & Management", status = "primary", solidHeader = TRUE,
            width = 12, height = "500px",
            
            conditionalPanel(
              condition = "!output.fileUploaded",
              div(class = "text-center", style = "padding: 80px;",
                  icon("upload", class = "fa-5x", style = "color: #ddd;"),
                  h2("Upload Expression Data", style = "color: #999; margin-top: 30px;"),
                  p("Drag & drop or browse for your gene expression matrix", 
                    style = "color: #999; font-size: 1.2em;"),
                  p("Supported: CSV, TXT, TSV, Excel", 
                    style = "color: #999; font-style: italic;")
              )
            ),
            
            conditionalPanel(
              condition = "output.fileUploaded",
              tabsetPanel(
                tabPanel("📋 Data Preview",
                  withSpinner(
                    DT::dataTableOutput("dataPreview"),
                    type = 1, color = "#3c8dbc"
                  )
                ),
                
                tabPanel("📈 Data Distribution",
                  fluidRow(
                    column(6, withSpinner(plotlyOutput("expressionDist"), type = 1)),
                    column(6, withSpinner(plotlyOutput("sampleCorrelation"), type = 1))
                  )
                ),
                
                tabPanel("🔍 Missing Data",
                  withSpinner(plotOutput("missingDataPlot"), type = 1)
                )
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.fileUploaded",
          fluidRow(
            valueBoxOutput("nGenes", width = 3),
            valueBoxOutput("nSamples", width = 3), 
            valueBoxOutput("dataSize", width = 3),
            valueBoxOutput("dataCompleteness", width = 3)
          ),
          
          fluidRow(
            box(
              title = "📊 Data Quality Summary", status = "info", solidHeader = TRUE,
              width = 12,
              div(class = "status-card",
                  fluidRow(
                    column(4, 
                      h4("Expression Range", style = "color: #3c8dbc;"),
                      textOutput("expressionRange")
                    ),
                    column(4,
                      h4("Sample Quality", style = "color: #3c8dbc;"),
                      textOutput("sampleQuality")
                    ),
                    column(4,
                      h4("Normalization Status", style = "color: #3c8dbc;"),
                      textOutput("normalizationStatus")
                    )
                  )
              )
            )
          )
        )
      ),
      
      # Data Filtering Tab
      tabItem(tabName = "filtering",
        fluidRow(
          box(
            title = "🔍 Advanced Data Filtering", status = "warning", solidHeader = TRUE,
            width = 8,
            conditionalPanel(
              condition = "output.fileUploaded",
              tabsetPanel(
                tabPanel("🧬 Gene Filtering",
                  div(class = "custom-tab-content",
                      h4("Expression-based Filtering"),
                      fluidRow(
                        column(6,
                          sliderInput("minExpressionLevel", "Min Expression:",
                                     min = 0, max = 20, value = 1, step = 0.1)
                        ),
                        column(6,
                          sliderInput("minExpressingSamples", "Min Expressing Samples:",
                                     min = 1, max = 100, value = 10, step = 1)
                        )
                      ),
                      
                      h4("Variance-based Filtering"),
                      sliderInput("varianceThreshold", "Keep top % by variance:",
                                 min = 10, max = 100, value = 75, step = 5),
                      
                      h4("Advanced Filters"),
                      checkboxInput("removeConstant", "Remove constant genes", value = TRUE),
                      checkboxInput("removeLowVar", "Remove low-variance genes", value = TRUE)
                  )
                ),
                
                tabPanel("👥 Sample QC",
                  div(class = "custom-tab-content",
                      h4("Sample Quality Control"),
                      checkboxInput("enableSampleQC", "Enable automatic sample QC", value = TRUE),
                      
                      conditionalPanel(
                        condition = "input.enableSampleQC",
                        sliderInput("outlierSDThreshold", "Outlier SD threshold:",
                                   min = 1, max = 5, value = 2.5, step = 0.1),
                        
                        sliderInput("correlationThreshold", "Min sample correlation:",
                                   min = 0.5, max = 1, value = 0.8, step = 0.05),
                        
                        checkboxInput("showOutliers", "Highlight outlier samples", value = TRUE)
                      )
                  )
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.fileUploaded",
              div(class = "text-center", style = "padding: 100px;",
                  icon("filter", class = "fa-4x", style = "color: #ccc;"),
                  h3("Upload data to enable filtering", style = "color: #999;")
              )
            )
          ),
          
          box(
            title = "📊 Filter Results", status = "success", solidHeader = TRUE,
            width = 4,
            conditionalPanel(
              condition = "output.fileUploaded",
              div(class = "status-card",
                  h4("Before Filtering", style = "color: #3c8dbc;"),
                  p(paste("Genes:", textOutput("originalGenes", inline = TRUE))),
                  p(paste("Samples:", textOutput("originalSamples", inline = TRUE))),
                  
                  hr(),
                  h4("After Filtering", style = "color: #27ae60;"),
                  p(paste("Genes:", textOutput("filteredGenes", inline = TRUE))),
                  p(paste("Samples:", textOutput("filteredSamples", inline = TRUE))),
                  
                  hr(),
                  actionBttn("applyFilters", "🔄 Apply Filters",
                            style = "gradient", color = "primary", size = "lg", block = TRUE)
              )
            )
          )
        )
      ),

      # Sample QC Tab  
      tabItem(tabName = "sample_qc",
        fluidRow(
          box(
            title = "👥 Sample Quality Control Dashboard", status = "info", solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.fileUploaded",
              tabsetPanel(
                tabPanel("📊 Sample Clustering",
                  withSpinner(plotlyOutput("sampleClusterPlot", height = "500px"), type = 1)
                ),
                
                tabPanel("🎯 Outlier Detection", 
                  fluidRow(
                    column(6, withSpinner(plotlyOutput("outlierScatterPlot"), type = 1)),
                    column(6, withSpinner(plotOutput("outlierBoxPlot"), type = 1))
                  )
                ),
                
                tabPanel("🔗 Sample Correlations",
                  withSpinner(plotlyOutput("sampleCorHeatmap", height = "600px"), type = 1)
                )
              )
            )
          )
        )
      ),

      # Network Analysis Tab Enhanced
      tabItem(tabName = "network",
        fluidRow(
          box(
            title = "🌐 Network Construction & Topology Analysis", status = "primary", solidHeader = TRUE,
            width = 12, height = "700px",
            
            conditionalPanel(
              condition = "!output.analysisComplete",
              div(class = "text-center", style = "padding: 150px;",
                  icon("network-wired", class = "fa-5x", style = "color: #ddd;"),
                  h2("Execute Analysis Pipeline", style = "color: #999; margin-top: 30px;"),
                  p("Configure parameters and run WGCNA to explore network topology", 
                    style = "color: #999; font-size: 1.2em;")
              )
            ),
            
            conditionalPanel(
              condition = "output.analysisComplete",
              tabsetPanel(
                tabPanel("⚡ Power Selection",
                  fluidRow(
                    column(8, withSpinner(plotlyOutput("powerSelectionPlot", height = "500px"), type = 1)),
                    column(4,
                      div(class = "status-card",
                          h4("Selected Parameters", style = "color: #3c8dbc;"),
                          verbatimTextOutput("powerSummary"),
                          
                          hr(),
                          h4("Network Quality", style = "color: #3c8dbc;"),
                          div(class = "metric-box", style = "margin: 10px 0;",
                              div(class = "metric-value", textOutput("networkR2")),
                              div(class = "metric-label", "Scale-Free R²")
                          )
                      )
                    )
                  )
                ),
                
                tabPanel("📊 Connectivity Analysis",
                  fluidRow(
                    column(6, withSpinner(plotlyOutput("connectivityPlot"), type = 1)),
                    column(6, withSpinner(plotlyOutput("scaleFreePlot"), type = 1))
                  ),
                  fluidRow(
                    column(12, withSpinner(plotlyOutput("meanConnectivityPlot"), type = 1))
                  )
                ),
                
                tabPanel("🔥 TOM Heatmap",
                  withSpinner(plotlyOutput("tomHeatmap", height = "600px"), type = 1)
                ),
                
                tabPanel("🌐 Network Graph",
                  withSpinner(plotlyOutput("networkGraphPlot", height = "600px"), type = 1)
                )
              )
            )
          )
        )
      ),

      # Enhanced Module Analysis Tab
      tabItem(tabName = "modules",
        fluidRow(
          box(
            title = "🎨 Module Detection & Characterization", status = "success", solidHeader = TRUE,
            width = 12, height = "750px",
            
            conditionalPanel(
              condition = "!output.analysisComplete",
              div(class = "text-center", style = "padding: 150px;",
                  icon("sitemap", class = "fa-5x", style = "color: #ddd;"),
                  h2("Module Analysis Pending", style = "color: #999; margin-top: 30px;"),
                  p("Complete network analysis to explore gene modules", 
                    style = "color: #999; font-size: 1.2em;")
              )
            ),
            
            conditionalPanel(
              condition = "output.analysisComplete",
              tabsetPanel(
                tabPanel("🌳 Dendrogram",
                  withSpinner(plotlyOutput("dendrogramPlot", height = "650px"), type = 1)
                ),
                
                tabPanel("🎨 Module Colors",
                  fluidRow(
                    column(8, withSpinner(plotOutput("moduleColorPlot", height = "500px"), type = 1)),
                    column(4,
                      div(class = "status-card",
                          h4("Module Statistics", style = "color: #3c8dbc;"),
                          textOutput("moduleStats"),
                          
                          hr(),
                          h4("Largest Modules", style = "color: #3c8dbc;"),
                          tableOutput("topModulesTable")
                      )
                    )
                  )
                ),
                
                tabPanel("🔗 Module Relationships",
                  withSpinner(plotlyOutput("moduleEigengenePlot", height = "600px"), type = 1)
                ),
                
                tabPanel("📊 Module Networks",
                  fluidRow(
                    column(6,
                      selectInput("selectedModule", "Select Module:",
                                 choices = NULL, selected = NULL)
                    ),
                    column(6,
                      numericInput("topGenes", "Top connected genes:", 
                                  value = 20, min = 5, max = 100, step = 5)
                    )
                  ),
                  withSpinner(plotlyOutput("moduleNetworkPlot", height = "500px"), type = 1)
                )
              )
            )
        )
      ),
      
      conditionalPanel(
        condition = "output.analysisComplete",
        fluidRow(
          box(
            title = "📋 Module Summary Table", status = "info", solidHeader = TRUE,
            width = 12,
            withSpinner(DT::dataTableOutput("moduleSummaryTable"), type = 1)
          )
        )
      )
    ),

    # New Gene Explorer Tab
    tabItem(tabName = "genes",
      fluidRow(
        box(
          title = "🔍 Gene-Level Analysis", status = "warning", solidHeader = TRUE,
          width = 8,
          conditionalPanel(
            condition = "output.analysisComplete",
            tabsetPanel(
              tabPanel("🎯 Gene Search",
                fluidRow(
                  column(6,
                    selectizeInput("searchGenes", "Search Genes:",
                                  choices = NULL, multiple = TRUE,
                                  options = list(
                                    placeholder = "Type gene names...",
                                    maxItems = 10
                                  ))
                  ),
                  column(6,
                    selectInput("searchModule", "Filter by Module:",
                               choices = NULL)
                  )
                ),
                
                conditionalPanel(
                  condition = "input.searchGenes.length > 0",
                  withSpinner(DT::dataTableOutput("geneInfoTable"), type = 1)
                )
              ),
              
              tabPanel("📈 Expression Profiles",
                conditionalPanel(
                  condition = "input.searchGenes.length > 0",
                  withSpinner(plotlyOutput("geneExpressionPlot", height = "500px"), type = 1)
                )
              ),
              
              tabPanel("🔗 Co-expression Network",
                conditionalPanel(
                  condition = "input.searchGenes.length > 0",
                  withSpinner(plotlyOutput("geneNetworkPlot", height = "500px"), type = 1)
                )
              )
            )
          ),
          
          conditionalPanel(
            condition = "!output.analysisComplete",
            div(class = "text-center", style = "padding: 100px;",
                icon("search", class = "fa-4x", style = "color: #ccc;"),
                h3("Complete analysis to explore genes", style = "color: #999;")
            )
          )
        ),
        
        box(
          title = "📊 Gene Statistics", status = "primary", solidHeader = TRUE,
          width = 4,
          conditionalPanel(
            condition = "output.analysisComplete",
            div(class = "status-card",
                h4("Module Distribution", style = "color: #3c8dbc;"),
                withSpinner(plotOutput("moduleDistributionPlot", height = "200px"), type = 1),
                
                hr(),
                h4("Top Hub Genes", style = "color: #3c8dbc;"),
                tableOutput("hubGenesTable")
            )
          )
        )
      )
    ),

    # Enhanced Quality Control Tab
    tabItem(tabName = "qc",
      fluidRow(
        box(
          title = "📈 Comprehensive Quality Metrics", status = "danger", solidHeader = TRUE,
          width = 12,
          conditionalPanel(
            condition = "output.analysisComplete",
            
            # Quality metrics overview
            fluidRow(
              column(3,
                div(class = "metric-box",
                    div(class = "metric-value", textOutput("selectedPower")),
                    div(class = "metric-label", "Selected Power")
                )
              ),
              column(3,
                div(class = "metric-box",
                    div(class = "metric-value", textOutput("scaleFreeFit")),
                    div(class = "metric-label", "Scale-Free R²")
                )
              ),
              column(3,
                div(class = "metric-box",
                    div(class = "metric-value", textOutput("meanConnectivity")),
                    div(class = "metric-label", "Mean Connectivity")
                )
              ),
              column(3,
                div(class = "metric-box",
                    div(class = "metric-value", textOutput("numModules")),
                    div(class = "metric-label", "Modules Detected")
                )
              )
            ),
            
            br(),
            tabsetPanel(
              tabPanel("📊 Network Quality",
                fluidRow(
                  column(6, withSpinner(plotlyOutput("networkQualityPlot"), type = 1)),
                  column(6, withSpinner(plotlyOutput("moduleStabilityPlot"), type = 1))
                )
              ),
              
              tabPanel("👥 Sample QC",
                fluidRow(
                  column(6, withSpinner(plotOutput("sampleQCPlot"), type = 1)),
                  column(6, withSpinner(plotOutput("sampleClusteringQC"), type = 1))
                )
              ),
              
              tabPanel("🧬 Gene QC", 
                fluidRow(
                  column(6, withSpinner(plotOutput("geneQCPlot"), type = 1)),
                  column(6, withSpinner(plotOutput("geneVariancePlot"), type = 1))
                )
              )
            )
          ),
          
          conditionalPanel(
            condition = "!output.analysisComplete",
            div(class = "text-center", style = "padding: 100px;",
                icon("chart-line", class = "fa-4x", style = "color: #ccc;"),
                h3("Execute analysis to view quality metrics", style = "color: #999;")
            )
          )
        )
      )
    ),

    # New Comparison Tab
    tabItem(tabName = "compare",
      fluidRow(
        box(
          title = "📊 Module & Analysis Comparisons", status = "info", solidHeader = TRUE,
          width = 12,
          conditionalPanel(
            condition = "output.analysisComplete",
            tabsetPanel(
              tabPanel("⚖️ Module Overlap",
                p("Compare modules across different parameter settings"),
                withSpinner(plotlyOutput("moduleOverlapPlot", height = "500px"), type = 1)
              ),
              
              tabPanel("📈 Parameter Sensitivity",
                p("Analyze how different parameters affect module detection"),
                withSpinner(plotlyOutput("parameterSensitivityPlot", height = "500px"), type = 1)
              ),
              
              tabPanel("🎯 Module Enrichment",
                p("Placeholder for gene set enrichment analysis"),
                div(class = "alert alert-info",
                    icon("info-circle"), " Gene set enrichment analysis coming soon!")
              )
            )
          ),
          
          conditionalPanel(
            condition = "!output.analysisComplete",
            div(class = "text-center", style = "padding: 100px;",
                icon("balance-scale", class = "fa-4x", style = "color: #ccc;"),
                h3("Complete analysis to enable comparisons", style = "color: #999;")
            )
          )
        )
      )
    ),

    # Enhanced Export Tab
    tabItem(tabName = "export",
      fluidRow(
        box(
          title = "📋 Advanced Export Center", status = "success", solidHeader = TRUE,
          width = 8,
          conditionalPanel(
            condition = "output.analysisComplete",
            
            h3("📁 Available Downloads"),
            tabsetPanel(
              tabPanel("📊 Data Exports",
                div(class = "custom-tab-content",
                    h4("Analysis Results"),
                    fluidRow(
                      column(4,
                        downloadBttn("downloadModules", "Module Assignments",
                                    style = "gradient", color = "primary", size = "md", block = TRUE)
                      ),
                      column(4,
                        downloadBttn("downloadNetwork", "Network Data",
                                    style = "gradient", color = "success", size = "md", block = TRUE)
                      ),
                      column(4,
                        downloadBttn("downloadExpression", "Filtered Expression",
                                    style = "gradient", color = "info", size = "md", block = TRUE)
                      )
                    ),
                    
                    br(),
                    h4("Connectivity Data"),
                    fluidRow(
                      column(6,
                        downloadBttn("downloadTOM", "TOM Matrix",
                                    style = "gradient", color = "warning", size = "md", block = TRUE)
                      ),
                      column(6,
                        downloadBttn("downloadEigengenes", "Module Eigengenes",
                                    style = "gradient", color = "danger", size = "md", block = TRUE)
                      )
                    )
                )
              ),
              
              tabPanel("📈 Visualization Exports",
                div(class = "custom-tab-content",
                    h4("Individual Plots"),
                    fluidRow(
                      column(4,
                        downloadBttn("downloadDendrogram", "Dendrogram",
                                    style = "gradient", color = "primary", size = "md", block = TRUE)
                      ),
                      column(4,
                        downloadBttn("downloadHeatmap", "TOM Heatmap", 
                                    style = "gradient", color = "success", size = "md", block = TRUE)
                      ),
                      column(4,
                        downloadBttn("downloadNetworkPlot", "Network Plot",
                                    style = "gradient", color = "info", size = "md", block = TRUE)
                      )
                    ),
                    
                    br(),
                    h4("Comprehensive Reports"),
                    fluidRow(
                      column(6,
                        downloadBttn("downloadAllPlots", "All Plots (PDF)",
                                    style = "gradient", color = "warning", size = "lg", block = TRUE)
                      ),
                      column(6,
                        downloadBttn("downloadReport", "Full Analysis Report",
                                    style = "gradient", color = "danger", size = "lg", block = TRUE)
                      )
                    )
                )
              ),
              
              tabPanel("⚙️ Analysis Settings",
                div(class = "custom-tab-content",
                    h4("Reproducibility"),
                    fluidRow(
                      column(6,
                        downloadBttn("downloadSettings", "Analysis Parameters",
                                    style = "gradient", color = "primary", size = "md", block = TRUE)
                      ),
                      column(6,
                        downloadBttn("downloadSessionInfo", "Session Info",
                                    style = "gradient", color = "success", size = "md", block = TRUE)
                      )
                    )
                )
              )
            )
          ),
          
          conditionalPanel(
            condition = "!output.analysisComplete",
            div(class = "text-center", style = "padding: 100px;",
                icon("download", class = "fa-4x", style = "color: #ccc;"),
                h2("Complete Analysis First", style = "color: #999; margin-top: 30px;"),
                p("Export options will be available after running WGCNA analysis", 
                  style = "color: #999; font-size: 1.2em;")
            )
          )
        ),
        
        box(
          title = "📊 Export Summary", status = "info", solidHeader = TRUE,
          width = 4,
          conditionalPanel(
            condition = "output.analysisComplete",
            div(class = "status-card",
                h4("Analysis Summary", style = "color: #3c8dbc;"),
                verbatimTextOutput("analysisSummary"),
                
                hr(),
                h4("Export Options", style = "color: #3c8dbc;"),
                checkboxInput("includeRawData", "Include raw data", value = TRUE),
                checkboxInput("includeParameters", "Include parameters", value = TRUE),
                checkboxInput("includeSessionInfo", "Include session info", value = TRUE),
                
                hr(),
                selectInput("exportFormat", "Default format:",
                           choices = list("CSV" = "csv", "Excel" = "xlsx", "RDS" = "rds"),
                           selected = "csv")
            )
          )
        )
      )
    ),

    # Enhanced Help Tab
    tabItem(tabName = "help",
      fluidRow(
        box(
          title = "📚 Documentation & User Guide", status = "info", solidHeader = TRUE,
          width = 12,
          
          tabsetPanel(
            tabPanel("🚀 Quick Start",
              div(class = "custom-tab-content",
                  h3("🧬 Welcome to Advanced WGCNA Explorer"),
                  p("This application provides a comprehensive interface for Weighted Gene Co-expression Network Analysis with advanced features and interactive visualizations."),
                  
                  h4("📋 Analysis Workflow:"),
                  tags$ol(
                    tags$li(strong("📊 Data Management:"), " Upload expression data and apply quality filters"),
                    tags$li(strong("⚙️ Parameter Configuration:"), " Set preprocessing and network construction parameters"),
                    tags$li(strong("🚀 Pipeline Execution:"), " Run the complete WGCNA analysis"),
                    tags$li(strong("🎨 Results Exploration:"), " Visualize modules and network topology"),
                    tags$li(strong("🔍 Gene Analysis:"), " Explore individual genes and their connections"),
                    tags$li(strong("📋 Export Results:"), " Download data, plots, and comprehensive reports")
                  )
              )
            ),
            
            tabPanel("📁 Data Requirements",
              div(class = "custom-tab-content",
                  h4("Expression Data Format:"),
                  tags$ul(
                    tags$li("Matrix format: Genes (rows) × Samples (columns)"),
                    tags$li("First column: Gene names/IDs (unique)"),
                    tags$li("Header row: Sample names (unique)"),
                    tags$li("Numeric expression values (no missing values in gene names)"),
                    tags$li("Recommended: Log-transformed or normalized expression")
                  ),
                  
                  h4("Supported File Formats:"),
                  tags$ul(
                    tags$li("CSV (Comma-separated values)"),
                    tags$li("TXT/TSV (Tab-separated values)"),
                    tags$li("Excel (.xlsx)")
                  ),
                  
                  h4("Sample Size Recommendations:"),
                  tags$ul(
                    tags$li("Minimum: 15-20 samples per group"),
                    tags$li("Recommended: 30+ samples for robust networks"),
                    tags$li("Genes: 5,000-20,000 for optimal performance")
                  )
              )
            ),
            
            tabPanel("⚙️ Parameters Guide",
              div(class = "custom-tab-content",
                  h4("🔧 Key Parameters:"),
                  
                  h5("Network Construction:"),
                  tags$ul(
                    tags$li(strong("Soft Power:"), " Controls network connectivity (higher = more stringent)"),
                    tags$li(strong("Network Type:"), " Unsigned (all correlations) vs Signed (direction matters)"),
                    tags$li(strong("Correlation Method:"), " Pearson vs Biweight midcorrelation (robust)")
                  ),
                  
                  h5("Module Detection:"),
                  tags$ul(
                    tags$li(strong("Deep Split:"), " Tree cutting sensitivity (0-4, higher = more modules)"),
                    tags$li(strong("Min Module Size:"), " Minimum genes per module (20-100 typical)"),
                    tags$li(strong("Merge Height:"), " Correlation threshold for merging similar modules")
                  ),
                  
                  h5("Quality Thresholds:"),
                  tags$ul(
                    tags$li(strong("Scale-Free R²:"), " Network topology fit (>0.8 recommended)"),
                    tags$li(strong("Mean Connectivity:"), " Average gene connections (<200 typical)")
                  )
              )
            ),
            
            tabPanel("🔗 Resources",
              div(class = "custom-tab-content",
                  h4("📖 Documentation:"),
                  tags$ul(
                    tags$li(tags$a("WGCNA Official Tutorial", 
                                  href = "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/", 
                                  target = "_blank")),
                    tags$li(tags$a("WGCNA R Package", 
                                  href = "https://cran.r-project.org/package=WGCNA", 
                                  target = "_blank")),
                    tags$li(tags$a("Co-expression Network Analysis Principles", 
                                  href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559", 
                                  target = "_blank"))
                  ),
                  
                  h4("🎓 Learning Resources:"),
                  tags$ul(
                    tags$li(tags$a("WGCNA Workshop Materials", href = "#", target = "_blank")),
                    tags$li(tags$a("Gene Co-expression Analysis Guide", href = "#", target = "_blank")),
                    tags$li(tags$a("Network Biology Fundamentals", href = "#", target = "_blank"))
                  ),
                  
                  h4("💬 Support:"),
                  tags$ul(
                    tags$li(tags$a("GitHub Issues", href = "https://github.com/yourusername/MyWGCNAResearchProject/issues", target = "_blank")),
                    tags$li(tags$a("WGCNA Google Group", href = "https://groups.google.com/g/wgcna", target = "_blank"))
                  )
              )
            )
          )
        )
      )
    )
  )
)