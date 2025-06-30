# Enhanced WGCNA Explorer UI
# ============================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)

# Define Enhanced UI
ui <- dashboardPage(
  skin = "blue",
  
  # Enhanced Header
  dashboardHeader(
    title = "🧬 Enhanced WGCNA Explorer",
    titleWidth = 350
  ),
  
  # Enhanced Sidebar
  dashboardSidebar(
    width = 280,
    
    sidebarMenu(
      id = "tabs",
      menuItem("📊 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("⚙️ Settings", tabName = "settings", icon = icon("cogs")),
      menuItem("🌳 Network", tabName = "network", icon = icon("project-diagram")),
      menuItem("🎨 Modules", tabName = "modules", icon = icon("palette")),
      menuItem("🔍 Genes", tabName = "genes", icon = icon("search")),
      menuItem("📋 Export", tabName = "export", icon = icon("download")),
      menuItem("ℹ️ Help", tabName = "help", icon = icon("info-circle"))
    )
  ),

  # Enhanced Dashboard Body
  dashboardBody(
    useShinyjs(),
    
    # Enhanced CSS styling
    tags$head(
      tags$style(HTML("
        .main-header .navbar {
          background-color: #3c8dbc;
        }
        .main-header .logo {
          background-color: #367fa9;
        }
        .box {
          border-radius: 5px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .btn-primary {
          background-color: #3c8dbc;
          border-color: #367fa9;
        }
        .btn-primary:hover {
          background-color: #367fa9;
          border-color: #2e6da4;
        }
      "))
    ),

    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "📊 Data Management", status = "primary", solidHeader = TRUE,
            width = 12,
            
            fileInput("file", "Choose Expression Data File",
                     accept = c(".csv", ".txt", ".xlsx", ".xls")),
            
            conditionalPanel(
              condition = "output.fileUploaded",
              DT::dataTableOutput("dataTable") %>% withSpinner(color = "#3c8dbc")
            ),
            
            conditionalPanel(
              condition = "!output.fileUploaded",
              div(class = "text-center", style = "padding: 50px;",
                  h4("Upload your gene expression data to begin"),
                  p("Supported formats: CSV, TXT, Excel")
              )
            )
          )
        )
      ),
      
      # Settings Tab
      tabItem(tabName = "settings",
        fluidRow(
          box(
            title = "⚙️ Analysis Settings", status = "primary", solidHeader = TRUE,
            width = 6,
            
            h4("Network Parameters"),
            sliderInput("softPower", "Soft Threshold Power:",
                       min = 1, max = 20, value = 6, step = 1),
            
            selectInput("networkType", "Network Type:",
                       choices = list("unsigned" = "unsigned",
                                    "signed" = "signed",
                                    "signed hybrid" = "signed hybrid"),
                       selected = "signed"),
            
            h4("Module Detection"),
            sliderInput("minModuleSize", "Minimum Module Size:",
                       min = 10, max = 100, value = 30, step = 5),
            
            sliderInput("mergeCutHeight", "Module Merge Cut Height:",
                       min = 0.1, max = 0.5, value = 0.25, step = 0.05),
            
            br(),
            actionButton("runAnalysis", "Run WGCNA Analysis", 
                        class = "btn btn-success", style = "width: 100%;")
          ),
          
          box(
            title = "📈 Analysis Progress", status = "info", solidHeader = TRUE,
            width = 6,
            
            conditionalPanel(
              condition = "output.analysisRunning",
              div(style = "text-align: center; padding: 30px;",
                  icon("spinner", class = "fa-spin", style = "font-size: 48px; color: #3c8dbc;"),
                  h4("Analysis in Progress..."),
                  textOutput("progressText")
              )
            ),
            
            conditionalPanel(
              condition = "!output.analysisRunning",
              div(style = "text-align: center; padding: 50px;",
                  h4("Ready to Run Analysis"),
                  p("Configure parameters and click 'Run WGCNA Analysis'")
              )
            )
          )
        )
      ),
      
      # Network Tab
      tabItem(tabName = "network",
        fluidRow(
          box(
            title = "🌳 Network Visualization", status = "primary", solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.networkAvailable",
              tabsetPanel(
                tabPanel("Network Heatmap",
                         plotlyOutput("networkHeatmap") %>% withSpinner(color = "#3c8dbc")
                ),
                tabPanel("Connectivity",
                         plotlyOutput("connectivityPlot") %>% withSpinner(color = "#3c8dbc")
                ),
                tabPanel("Network Stats",
                         verbatimTextOutput("networkStats")
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.networkAvailable",
              div(class = "text-center", style = "padding: 100px;",
                  h4("Run Analysis to View Network"),
                  p("Complete the WGCNA analysis to explore network properties.")
              )
            )
          )
        )
      ),
      
      # Modules Tab
      tabItem(tabName = "modules",
        fluidRow(
          box(
            title = "🎨 Module Analysis", status = "primary", solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.modulesAvailable",
              tabsetPanel(
                tabPanel("Module Dendrogram",
                         plotlyOutput("moduleDendrogram") %>% withSpinner(color = "#3c8dbc")
                ),
                tabPanel("Module Colors",
                         plotlyOutput("moduleColors") %>% withSpinner(color = "#3c8dbc")
                ),
                tabPanel("Module Summary",
                         DT::dataTableOutput("moduleSummary") %>% withSpinner(color = "#3c8dbc")
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.modulesAvailable",
              div(class = "text-center", style = "padding: 100px;",
                  h4("Modules Not Yet Detected"),
                  p("Run WGCNA analysis to detect and visualize modules.")
              )
            )
          )
        )
      ),
      
      # Genes Tab
      tabItem(tabName = "genes",
        fluidRow(
          box(
            title = "🔍 Gene Explorer", status = "primary", solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.modulesAvailable",
              fluidRow(
                column(4,
                       selectizeInput("selectedGenes", "Search Genes:",
                                     choices = NULL, multiple = TRUE),
                       selectInput("selectedModule", "Select Module:",
                                  choices = NULL),
                       actionButton("analyzeGenes", "Analyze Selection",
                                   class = "btn btn-primary")
                ),
                column(8,
                       conditionalPanel(
                         condition = "output.geneAnalysisAvailable",
                         tabsetPanel(
                           tabPanel("Expression Plot",
                                    plotlyOutput("geneExpressionPlot") %>% withSpinner(color = "#3c8dbc")
                           ),
                           tabPanel("Module Membership",
                                    plotlyOutput("moduleMembershipPlot") %>% withSpinner(color = "#3c8dbc")
                           )
                         )
                       ),
                       conditionalPanel(
                         condition = "!output.geneAnalysisAvailable",
                         div(style = "text-align: center; padding: 80px;",
                             h4("Select genes to analyze"),
                             p("Choose genes or modules to explore.")
                         )
                       )
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.modulesAvailable",
              div(class = "text-center", style = "padding: 100px;",
                  h4("Complete Analysis First"),
                  p("Run WGCNA analysis to explore individual genes.")
              )
            )
          )
        )
      ),
      
      # Export Tab
      tabItem(tabName = "export",
        fluidRow(
          box(
            title = "📋 Export Results", status = "primary", solidHeader = TRUE,
            width = 12,
            
            conditionalPanel(
              condition = "output.resultsAvailable",
              fluidRow(
                column(4,
                       h4("Export Options"),
                       checkboxGroupInput("exportItems", "Select items to export:",
                                         choices = list(
                                           "Module assignments" = "modules",
                                           "Module eigengenes" = "eigengenes",
                                           "Network adjacency" = "adjacency",
                                           "Hub genes" = "hub_genes"
                                         ),
                                         selected = c("modules", "eigengenes")),
                       
                       selectInput("exportFormat", "File format:",
                                  choices = list("CSV" = "csv", "Excel" = "xlsx"),
                                  selected = "csv"),
                       
                       br(),
                       downloadButton("downloadResults", "Download Results",
                                     class = "btn btn-success", style = "width: 100%;")
                ),
                column(8,
                       h4("Export Preview"),
                       DT::dataTableOutput("exportPreview") %>% withSpinner(color = "#3c8dbc")
                )
              )
            ),
            
            conditionalPanel(
              condition = "!output.resultsAvailable",
              div(class = "text-center", style = "padding: 100px;",
                  h4("No Results to Export"),
                  p("Complete WGCNA analysis to export results.")
              )
            )
          )
        )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
        fluidRow(
          box(
            title = "ℹ️ Documentation", status = "primary", solidHeader = TRUE,
            width = 12,
            
            h3("Enhanced WGCNA Explorer"),
            h4("🚀 Quick Start"),
            tags$ol(
              tags$li("Upload your gene expression data (CSV, TXT, or Excel format)"),
              tags$li("Configure analysis parameters in the Settings tab"),
              tags$li("Run WGCNA analysis and wait for completion"),
              tags$li("Explore network and modules in respective tabs"),
              tags$li("Export results when analysis is complete")
            ),
            
            br(),
            h4("📊 Data Requirements"),
            tags$ul(
              tags$li("Expression data with genes as rows, samples as columns"),
              tags$li("Minimum 15-20 samples recommended"),
              tags$li("Pre-processed, normalized expression values"),
              tags$li("No missing values or proper handling of missing data")
            ),
            
            br(),
            h4("⚙️ Parameter Guidelines"),
            tags$ul(
              tags$li("Soft threshold power: Usually 6-20, check scale-free topology"),
              tags$li("Network type: 'signed' preserves correlation direction"),
              tags$li("Module size: 30+ genes for stable modules"),
              tags$li("Merge height: 0.25 is typically appropriate")
            ),
            
            br(),
            h4("🔗 Resources"),
            tags$ul(
              tags$li(tags$a("WGCNA Tutorial", 
                            href = "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/", 
                            target = "_blank")),
              tags$li(tags$a("WGCNA Package", 
                            href = "https://cran.r-project.org/web/packages/WGCNA/", 
                            target = "_blank"))
            )
          )
        )
      )
    )
  )
)
