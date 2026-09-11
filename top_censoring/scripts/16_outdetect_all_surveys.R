#!/usr/bin/env Rscript

# Outlier detection (outdetect-style, after Belotti/Mancini/Vecchi) for every
# microdata survey in 01-input/country. Flags observations whose log-welfare
# is too far from the median in units of a robust, MAD-based scale, which is
# insensitive to the very outliers being screened for (unlike mean/SD).
#
# robust_z = (log(welfare) - median(log(welfare))) / (1.4826 * MAD(log(welfare)))
# alpha = 3 is the conventional cutoff (~ Hampel identifier).
# Both tails are flagged: bottom_outlier (robust_z < -alpha) and
# top_outlier (robust_z > alpha).

suppressPackageStartupMessages(library(haven))

input_dir <- file.path("01-input", "country")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 3

infer_weight <- function(columns) {
  for (candidate in c("weight_p", "weight", "weight_h")) {
    if (candidate %in% columns) return(candidate)
  }
  NA_character_
}

to_ppp_day <- function(d) {
  denominator <- d$cpi2021 * d$icp2021 * 365
  if (all(is.finite(denominator)) && all(denominator > 0)) {
    return(d$welfare / denominator)
  }
  stop("Missing or invalid PPP conversion fields.")
}

files <- sort(list.files(input_dir, pattern = "\\.dta$", full.names = TRUE))
results <- list()

for (path in files) {
  d <- read_dta(path)
  columns <- names(d)
  if ("bins" %in% columns || !"welfare" %in% columns) next

  weight_column <- infer_weight(columns)
  if (is.na(weight_column)) next

  welfare <- to_ppp_day(d)
  weight <- as.numeric(d[[weight_column]])
  valid <- is.finite(welfare) & welfare > 0 & is.finite(weight) & weight > 0
  welfare <- welfare[valid]
  weight <- weight[valid]
  if (!length(welfare)) next

  log_welfare <- log(welfare)
  median_log <- median(log_welfare)
  mad_log <- median(abs(log_welfare - median_log))
  scale_log <- 1.4826 * mad_log
  robust_z <- (log_welfare - median_log) / scale_log

  bottom_outlier <- robust_z < -alpha
  top_outlier <- robust_z > alpha
  any_outlier <- bottom_outlier | top_outlier

  results[[length(results) + 1]] <- data.frame(
    file = basename(path),
    country_code = if ("countrycode" %in% columns) as.character(d$countrycode[1]) else sub("_.*", "", basename(path)),
    survey_year = if ("year" %in% columns) as.integer(d$year[1]) else NA_integer_,
    survey_name = if ("survname" %in% columns) as.character(d$survname[1]) else NA_character_,
    weight_column = weight_column,
    valid_records = length(welfare),
    median_log = median_log,
    mad_log = mad_log,
    alpha = alpha,
    n_bottom_outliers = sum(bottom_outlier),
    n_top_outliers = sum(top_outlier),
    n_any_outliers = sum(any_outlier),
    share_bottom_outliers_unweighted_pct = 100 * mean(bottom_outlier),
    share_top_outliers_unweighted_pct = 100 * mean(top_outlier),
    share_any_outliers_unweighted_pct = 100 * mean(any_outlier),
    weighted_population_bottom_outliers = sum(weight[bottom_outlier]),
    weighted_population_top_outliers = sum(weight[top_outlier]),
    share_bottom_outliers_weighted_pct = 100 * sum(weight[bottom_outlier]) / sum(weight),
    share_top_outliers_weighted_pct = 100 * sum(weight[top_outlier]) / sum(weight),
    share_any_outliers_weighted_pct = 100 * sum(weight[any_outlier]) / sum(weight),
    max_robust_z = max(robust_z),
    min_robust_z = min(robust_z),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, results)
summary <- summary[order(summary$share_any_outliers_weighted_pct, decreasing = TRUE), ]
write.csv(
  summary,
  file.path(output_dir, "outdetect_all_surveys_summary.csv"),
  row.names = FALSE
)

cat("Wrote", nrow(summary), "survey rows to outdetect_all_surveys_summary.csv\n")
print(summary[, c("country_code", "survey_year", "survey_name", "valid_records",
                   "n_bottom_outliers", "n_top_outliers",
                   "share_any_outliers_weighted_pct")])
