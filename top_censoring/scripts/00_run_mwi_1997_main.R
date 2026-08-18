# One-command entry point for the main Malawi 1997 top-censoring analysis.
# Run this file from the bindata_check repository root.

source(file.path("top_censoring", "scripts", "01_lis_top_diagnostics.R"))
source(file.path("top_censoring", "scripts", "02_mwi_1997_pareto.R"))
source(file.path("top_censoring", "scripts", "12_mwi_1997_eusilc_approaches.R"))

# Convert the main long-form output to an Excel-friendly comparison table.
output_dir <- file.path("top_censoring", "outputs")
comparison <- read.csv(
  file.path(output_dir, "mwi_1997_eusilc_approaches_comparison.csv"),
  stringsAsFactors = FALSE
)
comparison_table <- reshape(
  comparison[, c("approach", "indicator", "value")],
  idvar = "approach",
  timevar = "indicator",
  direction = "wide"
)
names(comparison_table) <- sub("^value\\.", "", names(comparison_table))
write.csv(
  comparison_table,
  file.path(output_dir, "mwi_1997_eusilc_approaches_table.csv"),
  row.names = FALSE
)

message("Main Malawi 1997 analysis is complete.")
message("Open: top_censoring/outputs/mwi_1997_eusilc_approaches_table.csv")
