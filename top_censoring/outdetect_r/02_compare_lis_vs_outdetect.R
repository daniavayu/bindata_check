#!/usr/bin/env Rscript

# Compare the LIS Tukey-fence screening (existing method) against the
# outdetect-EXACT screening (median + S-estimator via robustbase::Sn(),
# log space, alpha=3, top-tail only), over the same ~1,825 survey-year
# universe from GMD_all_2017.dta. Mirrors ../outdetect/02_compare_lis_vs_outdetect.py.
#
# Outputs:
#   - lis_vs_outdetect_comparison.csv : full universe, both methods side by
#     side, with flags for which method(s) singled out each survey.
#   - always_problematic_surveys.csv : the subset flagged by BOTH methods.

LIS_FULL <- file.path("top_censoring", "outputs", "all_surveys_four_indicators_lis_comparison.csv")
LIS_TOP24 <- file.path("top_censoring", "outputs", "24_surveys_four_indicators_lis_comparison.csv")
OUTDETECT_FULL <- file.path("top_censoring", "outdetect_r", "outputs", "outdetect_all_surveys_summary.csv")
WELFARE_TYPE <- file.path("corrected_method", "01-processing", "raw", "interpolated_means.csv")
OUT_DIR <- file.path("top_censoring", "outdetect_r", "outputs")

lis_full <- read.csv(LIS_FULL, stringsAsFactors = FALSE, check.names = FALSE)
lis_top24 <- read.csv(LIS_TOP24, stringsAsFactors = FALSE, check.names = FALSE)
outdetect_full <- read.csv(OUTDETECT_FULL, stringsAsFactors = FALSE, check.names = FALSE)

lis_top24_key <- paste(lis_top24$Country, lis_top24$Year, lis_top24$Survey, sep = "||")
outdetect_full$`Flagged by outdetect` <- suppressWarnings(as.numeric(outdetect_full[["Welfare share flagged (%)"]])) >= 5

merged <- merge(
  lis_full, outdetect_full,
  by = c("Country", "Year", "Survey"),
  suffixes = c(" (LIS)", " (outdetect)"),
  all = TRUE
)

merged$`Flagged by LIS` <- paste(merged$Country, merged$Year, merged$Survey, sep = "||") %in% lis_top24_key
merged$`Flagged by both (always problematic)` <- merged$`Flagged by LIS` & merged$`Flagged by outdetect`

welfare_type <- read.csv(WELFARE_TYPE, stringsAsFactors = FALSE)[, c("country_code", "year", "welfare_type")]
welfare_type <- welfare_type[!duplicated(welfare_type[, c("country_code", "year")]), ]
merged <- merge(
  merged, welfare_type,
  by.x = c("Country", "Year"), by.y = c("country_code", "year"),
  all.x = TRUE
)

lis_share_col <- grep("^Welfare share above ceiling", names(merged), value = TRUE)[1]
if (is.na(lis_share_col)) stop("Could not find the LIS welfare-share column after merge.")
lis_share <- suppressWarnings(as.numeric(merged[[lis_share_col]]))
merged <- merged[order(-merged$`Flagged by both (always problematic)`, -lis_share), ]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
comparison_path <- file.path(OUT_DIR, "lis_vs_outdetect_comparison.csv")
write.csv(merged, comparison_path, row.names = FALSE)

always_problematic <- merged[merged$`Flagged by both (always problematic)`, ]
always_problematic_path <- file.path(OUT_DIR, "always_problematic_surveys.csv")
write.csv(always_problematic, always_problematic_path, row.names = FALSE)

message("Wrote ", nrow(merged), " rows to ", comparison_path)
message("Wrote ", nrow(always_problematic), " rows to ", always_problematic_path)
message("")
message("Surveys flagged as problematic by BOTH methods:")
print(
  always_problematic[, c(
    "Country", "Year", "Survey", "welfare_type",
    "Welfare share above ceiling (%)", "Welfare share flagged (%)"
  )],
  row.names = FALSE
)
message("")
message("LIS-only count: ", sum(merged$`Flagged by LIS` & !merged$`Flagged by outdetect`))
message("outdetect-only count: ", sum(!merged$`Flagged by LIS` & merged$`Flagged by outdetect`))
