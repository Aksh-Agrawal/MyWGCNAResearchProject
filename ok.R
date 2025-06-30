# ok.R — development entrypoint

# 1. Set your project root
setwd("d:/r/New folder (3)/MyWGCNAResearchProject")

# 2. Load devtools so we can source the package code
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::load_all("Rpackage")    # sources everything under Rpackage/R/

# 3. Launch the Shiny app
shiny::runApp("Rpackage/inst/shiny", launch.browser = TRUE)
