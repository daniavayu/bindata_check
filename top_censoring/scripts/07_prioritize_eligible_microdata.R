#!/usr/bin/env Rscript

# Create the final priority list from the first-pass screening. Files with fewer
# than 100 valid records are treated as grouped distributions, not microdata.

input_file <- file.path(
  "top_censoring", "outputs", "all_microdata_lis_top_screening.csv"
)
output_file <- file.path(
  "top_censoring", "outputs", "eligible_microdata_lis_priority.csv"
)

if (!file.exists(input_file)) {
  stop("Run 06_screen_all_microdata.R before running this script.")
}

screening <- read.csv(input_file, stringsAsFactors = FALSE)
eligible <- screening[screening$valid_records >= 100, ]
eligible <- eligible[order(
  eligible$welfare_share_above_lis_pct,
  decreasing = TRUE
), ]

write.csv(eligible, output_file, row.names = FALSE)
print(eligible, row.names = FALSE)
message("Wrote eligible microdata priority list to ", output_file)
