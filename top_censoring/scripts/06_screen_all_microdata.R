#!/usr/bin/env Rscript

# First-pass top-tail screening for every microdata survey in 01-input/country.
# Grouped distributions (files with a `bins` column) are excluded because they
# cannot identify individual extreme observations.

suppressPackageStartupMessages(library(haven))

input_dir <- file.path("01-input", "country")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

weighted_quantile <- function(x, w, p) {
  order_index <- order(x)
  sorted_x <- x[order_index]
  sorted_w <- w[order_index]
  sorted_x[which(cumsum(sorted_w) / sum(sorted_w) >= p)[1]]
}

weighted_mean <- function(x, w) sum(x * w) / sum(w)

weighted_gini <- function(x, w) {
  order_index <- order(x)
  x <- x[order_index]
  w <- w[order_index]
  population <- c(0, cumsum(w) / sum(w))
  welfare <- c(0, cumsum(x * w) / sum(x * w))
  1 - sum(diff(population) * (head(welfare, -1) + tail(welfare, -1)))
}

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

  q1_log <- weighted_quantile(log(welfare), weight, 0.25)
  q3_log <- weighted_quantile(log(welfare), weight, 0.75)
  lis_ceiling <- exp(q3_log + 3 * (q3_log - q1_log))
  above_lis <- welfare > lis_ceiling
  top_coded <- pmin(welfare, lis_ceiling)

  mean_raw <- weighted_mean(welfare, weight)
  mean_coded <- weighted_mean(top_coded, weight)
  gini_raw <- weighted_gini(welfare, weight)
  gini_coded <- weighted_gini(top_coded, weight)

  results[[length(results) + 1]] <- data.frame(
    file = basename(path),
    country_code = if ("countrycode" %in% columns) as.character(d$countrycode[1]) else sub("_.*", "", basename(path)),
    survey_year = if ("year" %in% columns) as.integer(d$year[1]) else NA_integer_,
    survey_name = if ("survname" %in% columns) as.character(d$survname[1]) else NA_character_,
    weight_column = weight_column,
    valid_records = length(welfare),
    lis_ceiling_ppp_day = lis_ceiling,
    records_above_lis = sum(above_lis),
    population_share_above_lis_pct = 100 * sum(weight[above_lis]) / sum(weight),
    welfare_share_above_lis_pct = 100 * sum(welfare[above_lis] * weight[above_lis]) / sum(welfare * weight),
    maximum_to_lis_ceiling = max(welfare) / lis_ceiling,
    mean_change_from_lis_top_coding_pct = 100 * (mean_coded / mean_raw - 1),
    gini_change_from_lis_top_coding = gini_coded - gini_raw,
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, results)
summary <- summary[order(summary$welfare_share_above_lis_pct, decreasing = TRUE), ]
write.csv(
  summary,
  file.path(output_dir, "all_microdata_lis_top_screening.csv"),
  row.names = FALSE
)

print(summary, row.names = FALSE)
message("Wrote screening results to ", output_dir)
