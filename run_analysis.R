# Run from this folder with the cleaned data frame `analysis_data` in memory.
# If the original complete analysis has already been run, first assign:
# analysis_data <- Final_results$data
if(!exists("analysis_data", inherits=FALSE))
  stop("Create analysis_data from the already cleaned NHANES data before running this file.")
source("NHANES_manuscript_analysis.R", encoding="UTF-8")
Manuscript_results <- run_nhanes_manuscript(analysis_data)
Main_tables <- Manuscript_results$publication$Main
Supplement_tables <- Manuscript_results$publication$Supplementary
Figure_index <- Manuscript_results$figure_index
list2env(Manuscript_results$plots, envir=environment())
print(Figure_index, row.names=FALSE)
# Figures are not drawn or saved until requested in call_results.R.
