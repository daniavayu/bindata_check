#!/usr/bin/env Rscript

# Executive summary (EXACT outdetect): before/after comparison for the
# surveys flagged by both methods, mirroring
# ../outdetect/03_executive_summary.py, but with the S-estimator computed
# EXACTLY everywhere (robustbase::Sn(), Croux-Rousseeuw 1992 O(n log n)
# algorithm) -- no subsampling, no RNG, anywhere in this script, including
# the GMD-wide Table 1 (unlike the Python version, which still subsamples
# Table 1 because it is a plain O(n^2) implementation).
#
# The Malawi 1997 flagship deep dive (Tables 2 and 3) is recalculated from
# country/MWI_1997.dta, same as the Python version, so it shares one source
# with the other Malawi-specific analyses.
#
# Outputs:
#   outputs/executive_summary_always_problematic.csv
#   outputs/executive_summary_mwi_1997_flagship_top_values.csv
#   outputs/mwi_1997_treatment_comparison_5_scenarios.csv

suppressPackageStartupMessages({
  library(haven)
  library(robustbase)
})

DATA_PATH <- file.path("01-input", "GMD_all_2017.dta")
COUNTRY_FLAGSHIP_PATH <- file.path("01-input", "country", "MWI_1997.dta")
OUT_DIR <- file.path("top_censoring", "outdetect_r", "outputs")
ALWAYS_PROBLEMATIC_PATH <- file.path(OUT_DIR, "always_problematic_surveys.csv")

ALPHA <- 3
S_SCALE_FACTOR <- 1.1926

s_estimator_exact <- function(values) as.numeric(Sn(values, constant = S_SCALE_FACTOR))

weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

weighted_mean <- function(x, w) sum(x * w) / sum(w)

weighted_gini <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  population <- c(0, cumsum(w) / sum(w))
  welfare <- c(0, cumsum(x * w) / sum(x * w))
  1 - sum(diff(population) * (head(welfare, -1) + tail(welfare, -1)))
}

indicators <- function(x, w) {
  m <- weighted_mean(x, w)
  list(
    `Mean (PPP/day)` = m,
    Gini = weighted_gini(x, w),
    MLD = log(m) - weighted_mean(log(x), w),
    `Atkinson(2)` = 1 - (1 / weighted_mean(1 / x, w)) / m
  )
}

always_problematic <- read.csv(ALWAYS_PROBLEMATIC_PATH, stringsAsFactors = FALSE)
wanted <- always_problematic[, c("Country", "Year", "Survey")]
wanted_key <- paste(wanted$Country, wanted$Year, wanted$Survey, sep = "||")

# --- Load only the microdata for the target surveys -----------------------
message("Reading ", DATA_PATH, " to extract the target surveys...")
d <- read_dta(
  DATA_PATH,
  col_select = c(welfare_dec, weight, code, year, survname, use_microdata, use_bin, use_groupdata)
)
country <- toupper(as.character(d$code))
survey <- as.character(d$survname)
year <- as.integer(d$year)
welfare <- as.numeric(d$welfare_dec)
weight <- as.numeric(d$weight)
use_microdata <- as.numeric(d$use_microdata)
use_bin <- as.numeric(d$use_bin)
use_groupdata <- as.numeric(d$use_groupdata)
rm(d)
gc()

valid <- !is.na(use_microdata) & use_microdata == 1 &
  !is.na(use_bin) & use_bin == 0 &
  !is.na(use_groupdata) & use_groupdata == 0 &
  is.finite(welfare) & welfare > 0 &
  is.finite(weight) & weight > 0
key <- paste(country, year, survey, sep = "||")
keep <- valid & key %in% wanted_key

microdata <- data.frame(
  Country = country[keep], Year = year[keep], Survey = survey[keep],
  welfare = welfare[keep], weight = weight[keep],
  stringsAsFactors = FALSE
)
rm(country, survey, year, welfare, weight, use_microdata, use_bin, use_groupdata, valid, key, keep)
gc()

# --- Table 1: wide comparison for the "always problematic" surveys --------
rows <- list()
for (i in seq_len(nrow(wanted))) {
  ctry <- wanted$Country[i]
  yr <- wanted$Year[i]
  srv <- wanted$Survey[i]
  data <- microdata[microdata$Country == ctry & microdata$Year == yr & microdata$Survey == srv, ]
  values <- data$welfare
  weights <- data$weight
  log_values <- log(values)

  q <- weighted_quantile(log_values, weights, c(0.25, 0.75))
  lis_ceiling <- exp(q[2] + 3 * (q[2] - q[1]))
  lis_flag <- values > lis_ceiling
  lis_capped <- pmin(values, lis_ceiling)

  median_log <- median(log_values)
  s_log <- s_estimator_exact(log_values)
  outdetect_threshold <- exp(median_log + ALPHA * s_log)
  outdetect_flag <- values > outdetect_threshold
  outdetect_capped <- pmin(values, outdetect_threshold)

  raw <- indicators(values, weights)
  lis_ind <- indicators(lis_capped, weights)
  outdetect_ind <- indicators(outdetect_capped, weights)

  row <- list(
    Country = ctry, Year = yr, Survey = srv,
    `Valid records` = length(values),
    `LIS ceiling (PPP/day)` = lis_ceiling,
    `outdetect threshold (PPP/day)` = outdetect_threshold,
    `Population affected (%) - LIS` = 100 * sum(weights[lis_flag]) / sum(weights),
    `Population affected (%) - outdetect` = 100 * sum(weights[outdetect_flag]) / sum(weights),
    `Welfare affected (%) - LIS` = 100 * sum(values[lis_flag] * weights[lis_flag]) / sum(values * weights),
    `Welfare affected (%) - outdetect` = 100 * sum(values[outdetect_flag] * weights[outdetect_flag]) / sum(values * weights)
  )
  for (name in names(raw)) {
    row[[paste(name, "- Raw")]] <- raw[[name]]
    row[[paste(name, "- LIS-capped")]] <- lis_ind[[name]]
    row[[paste(name, "- Change (%) LIS")]] <- 100 * (lis_ind[[name]] / raw[[name]] - 1)
    row[[paste(name, "- outdetect-capped")]] <- outdetect_ind[[name]]
    row[[paste(name, "- Change (%) outdetect")]] <- 100 * (outdetect_ind[[name]] / raw[[name]] - 1)
  }
  rows[[i]] <- as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}
summary_df <- do.call(rbind, rows)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(OUT_DIR, "executive_summary_always_problematic.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
message("Wrote ", nrow(summary_df), " rows to ", summary_path)

# --- Table 2: Malawi 1997 flagship deep dive (top 30 largest values) ------
country_df <- read_dta(COUNTRY_FLAGSHIP_PATH)
cv <- is.finite(country_df$welfare) & country_df$welfare > 0 &
  is.finite(country_df$weight) & country_df$weight > 0 &
  is.finite(country_df$cpi2021) & country_df$cpi2021 > 0 &
  is.finite(country_df$icp2021) & country_df$icp2021 > 0

flagship_values <- as.numeric(
  country_df$welfare[cv] / (country_df$cpi2021[cv] * country_df$icp2021[cv] * 365)
)
flagship_weights <- as.numeric(country_df$weight[cv])
flagship_log <- log(flagship_values)

fq <- weighted_quantile(flagship_log, flagship_weights, c(0.25, 0.75))
lis_ceiling <- exp(fq[2] + 3 * (fq[2] - fq[1]))
flagship_s <- s_estimator_exact(flagship_log)
outdetect_threshold <- exp(median(flagship_log) + ALPHA * flagship_s)
lis_flag <- flagship_values > lis_ceiling
outdetect_flag <- flagship_values > outdetect_threshold

order_idx <- order(flagship_values, decreasing = TRUE)[1:30]
flagship_table <- data.frame(
  `Rank (largest first)` = seq_len(30),
  `Welfare (PPP/day)` = flagship_values[order_idx],
  `Sampling weight` = flagship_weights[order_idx],
  `Flagged by LIS (> ceiling)` = lis_flag[order_idx],
  `Flagged by outdetect (> threshold)` = outdetect_flag[order_idx],
  `Flagged by both` = lis_flag[order_idx] & outdetect_flag[order_idx],
  check.names = FALSE
)

flagship_path <- file.path(OUT_DIR, "executive_summary_mwi_1997_flagship_top_values.csv")
con <- file(flagship_path, "w")
writeLines(sprintf("# LIS ceiling (PPP/day): %.4f", lis_ceiling), con)
writeLines(sprintf("# outdetect threshold (PPP/day): %.4f", outdetect_threshold), con)
write.csv(flagship_table, con, row.names = FALSE)
close(con)
message("Wrote flagship top-30 table to ", flagship_path)
print(summary_df[, c("Country", "Year", "Survey", "Population affected (%) - LIS", "Population affected (%) - outdetect")], row.names = FALSE)

# --- Table 3: Malawi 1997 five-scenario treatment comparison --------------
two_highest <- order(flagship_values, decreasing = TRUE)[1:2]
keep_not_two <- rep(TRUE, length(flagship_values))
keep_not_two[two_highest] <- FALSE

scenarios <- list(
  Unadjusted = indicators(flagship_values, flagship_weights),
  `Remove 2 highest` = indicators(flagship_values[keep_not_two], flagship_weights[keep_not_two]),
  `Remove above LIS cap` = indicators(flagship_values[!lis_flag], flagship_weights[!lis_flag]),
  `Remove above outdetect cap` = indicators(flagship_values[!outdetect_flag], flagship_weights[!outdetect_flag]),
  `Cap at LIS cap` = indicators(pmin(flagship_values, lis_ceiling), flagship_weights),
  `Cap at outdetect cap` = indicators(pmin(flagship_values, outdetect_threshold), flagship_weights)
)
scenario_rows <- lapply(names(scenarios), function(name) {
  as.data.frame(c(list(Approach = name), scenarios[[name]]), check.names = FALSE, stringsAsFactors = FALSE)
})
scenario_table <- do.call(rbind, scenario_rows)

scenario_path <- file.path("top_censoring", "outdetect_r", "outputs", "mwi_1997_treatment_comparison_5_scenarios.csv")
con <- file(scenario_path, "w")
writeLines(sprintf("# LIS cap (PPP/day): %.4f", lis_ceiling), con)
writeLines(sprintf("# outdetect cap (PPP/day, exact S-estimator via robustbase::Sn(), no subsampling): %.4f", outdetect_threshold), con)
writeLines(sprintf("# Records above LIS cap: %d; records above outdetect cap: %d", sum(lis_flag), sum(outdetect_flag)), con)
write.csv(scenario_table, con, row.names = FALSE)
close(con)
message("Wrote 5-scenario treatment comparison to ", scenario_path)
print(scenario_table, row.names = FALSE, digits = 4)
