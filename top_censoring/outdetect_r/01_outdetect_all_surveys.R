#!/usr/bin/env Rscript

# outdetect (EXACT): robust top-tail outlier screening for every GMD survey.
#
# Same method as ../outdetect/01_outdetect_all_surveys.py:
#   outdetect welfare, norm(ln) zscore(median, s) alpha(3) out(top)
# (Belotti, Mancini and Vecchi, "Outlier Detection for Welfare Analysis",
# World Bank Policy Research Working Paper) -- center = unweighted median of
# log(welfare), scale = unweighted S-estimator (Rousseeuw & Croux, 1993) of
# log(welfare), flag values more than alpha=3 scales above the median.
#
# The only difference from the Python version: the S-estimator here is
# computed EXACTLY, with no random subsampling, using robustbase::Sn(), which
# implements the Croux-Rousseeuw (1992) O(n log n) algorithm -- the same fast
# algorithm Stata's `outdetect` uses internally in Mata. Benchmarked at ~0.5s
# for n = 1,000,000 and ~0.01s for n = 20,000, so every one of the ~1,825
# GMD surveys (even the largest, ~5 million records) gets a deterministic,
# reproducible threshold; no RANDOM_SEED, no S_ESTIMATOR_SUBSAMPLE_SIZE.
# See README.md for the correctness check against the naive O(n^2) formula.

suppressPackageStartupMessages({
  library(data.table)
  library(robustbase)
})

# NOTE: haven::read_dta() has no chunked-reading mode, so pointing it directly
# at the 3.6 GB GMD_all_2017.dta file (even with col_select) takes tens of
# minutes just to parse. 00_export_microdata.py reuses pandas' fast chunked
# Stata reader (same filter logic as ../outdetect/01_outdetect_all_surveys.py)
# to produce a small pre-filtered CSV in seconds; this script reads that CSV
# and does the actual EXACT S-estimator computation (robustbase::Sn()) here.
# Run `python top_censoring/outdetect_r/00_export_microdata.py` first if the
# CSV below does not exist yet.
MICRODATA_PATH <- file.path("top_censoring", "outdetect_r", "outputs", "_microdata_filtered.csv")
OUTPUT_PATH <- file.path("top_censoring", "outdetect_r", "outputs", "outdetect_all_surveys_summary.csv")

ALPHA <- 3
S_SCALE_FACTOR <- 1.1926 # Rousseeuw & Croux (1993) consistency constant for S

s_estimator_exact <- function(values) {
  as.numeric(Sn(values, constant = S_SCALE_FACTOR))
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
    mean_ppp_day = m,
    gini = weighted_gini(x, w),
    mld = log(m) - weighted_mean(log(x), w),
    atkinson_2 = 1 - (1 / weighted_mean(1 / x, w)) / m
  )
}

if (!file.exists(MICRODATA_PATH)) {
  stop(
    "Missing ", MICRODATA_PATH, " -- run `python top_censoring/outdetect_r/00_export_microdata.py` first."
  )
}

message("Reading ", MICRODATA_PATH, " (pre-filtered microdata export)...")
micro <- fread(MICRODATA_PATH, showProgress = FALSE)
setDF(micro)
micro$Country <- as.character(micro$Country)
micro$Survey <- as.character(micro$Survey)
micro$Year <- as.integer(micro$Year)

message("Screening ", nrow(micro), " valid records across all surveys...")
group_key <- interaction(micro$Survey, micro$Country, micro$Year, drop = TRUE)
groups <- split(seq_len(nrow(micro)), group_key)

results <- vector("list", length(groups))
for (i in seq_along(groups)) {
  idx <- groups[[i]]
  values <- micro$welfare[idx]
  weights <- micro$weight[idx]
  log_values <- log(values)

  median_log <- median(log_values) # unweighted, per Belotti/Mancini/Vecchi
  s_log <- s_estimator_exact(log_values) # unweighted, exact S-estimator (Rousseeuw & Croux, 1993)
  threshold <- if (s_log > 0) exp(median_log + ALPHA * s_log) else Inf

  top_flag <- values > threshold
  capped <- pmin(values, threshold)
  original <- indicators(values, weights)
  capped_ind <- indicators(capped, weights)
  welfare_share_flagged <- if (any(top_flag)) {
    100 * sum(values[top_flag] * weights[top_flag]) / sum(values * weights)
  } else {
    0
  }

  results[[i]] <- data.frame(
    Country = micro$Country[idx[1]],
    Year = micro$Year[idx[1]],
    Survey = micro$Survey[idx[1]],
    `Valid records` = length(values),
    `Median log welfare (unweighted)` = median_log,
    `S estimate, log space (unweighted, 1.1926 scaled)` = s_log,
    `outdetect threshold (PPP/day)` = threshold,
    `Records flagged (top)` = sum(top_flag),
    `Population share flagged (%)` = 100 * sum(weights[top_flag]) / sum(weights),
    `Welfare share flagged (%)` = welfare_share_flagged,
    `mean_ppp_day original` = original$mean_ppp_day,
    `mean_ppp_day outdetect-capped` = capped_ind$mean_ppp_day,
    `mean_ppp_day change (%)` = 100 * (capped_ind$mean_ppp_day / original$mean_ppp_day - 1),
    `gini original` = original$gini,
    `gini outdetect-capped` = capped_ind$gini,
    `gini change (%)` = 100 * (capped_ind$gini / original$gini - 1),
    `mld original` = original$mld,
    `mld outdetect-capped` = capped_ind$mld,
    `mld change (%)` = 100 * (capped_ind$mld / original$mld - 1),
    `atkinson_2 original` = original$atkinson_2,
    `atkinson_2 outdetect-capped` = capped_ind$atkinson_2,
    `atkinson_2 change (%)` = 100 * (capped_ind$atkinson_2 / original$atkinson_2 - 1),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, results)
summary_df <- summary_df[order(-summary_df[["Welfare share flagged (%)"]]), ]

dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, OUTPUT_PATH, row.names = FALSE)
message("Wrote ", nrow(summary_df), " rows to ", OUTPUT_PATH)
print(head(summary_df[, c("Country", "Year", "Survey", "Records flagged (top)", "Welfare share flagged (%)")], 30), row.names = FALSE)
