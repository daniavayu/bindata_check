#!/usr/bin/env Rscript

# Malawi 1997: comparison of approaches for high-income observations.
#
# This script follows the structure of Van Kerm (2007, EU-SILC): compare
# trimming, winsorizing, and a Pareto-tail model.  The Pareto shape is
# estimated with laeken::thetaPDC, which is robust to tail contamination.
# This is an approximation of the paper's OBRE estimator, not an exact
# reproduction of OBRE.  Survey weights are used in all estimands and in
# the Pareto estimation.

suppressPackageStartupMessages({
  library(haven)
  library(laeken)
})

input_file <- file.path("01-input", "country", "MWI_1997.dta")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ---------------------------------------------------------------
wq <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

wmean <- function(x, w) sum(x * w) / sum(w)

wgini <- function(x, w) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  pop <- c(0, cumsum(w) / sum(w))
  inc <- c(0, cumsum(x * w) / sum(x * w))
  1 - sum(diff(pop) * (head(inc, -1) + tail(inc, -1)))
}

income_share <- function(x, w, lower, upper) {
  q_lo <- wq(x, w, lower)
  q_hi <- wq(x, w, upper)
  include <- x >= q_lo & x <= q_hi
  sum(x[include] * w[include]) / sum(x * w)
}

indicators <- function(x, w, poverty_line) {
  mu <- wmean(x, w)
  poor <- x < poverty_line
  gap <- numeric(length(x))
  gap[poor] <- (poverty_line - x[poor]) / poverty_line
  ede05 <- wmean(sqrt(x), w)^2
  ede1 <- exp(wmean(log(x), w))
  ede2 <- 1 / wmean(1 / x, w)
  p10 <- wq(x, w, .10)
  p90 <- wq(x, w, .90)
  s20 <- income_share(x, w, 0, .20)
  s80 <- income_share(x, w, .80, 1)

  c(
    mean_ppp_day = mu,
    median_ppp_day = wq(x, w, .50),
    headcount_2_41 = sum(w[poor]) / sum(w),
    fgt1_2_41 = wmean(gap, w),
    fgt2_2_41 = wmean(gap^2, w),
    watts_2_41 = wmean(ifelse(poor, log(poverty_line / x), 0), w),
    p90_p10 = p90 / p10,
    s80_s20 = s80 / s20,
    gini = wgini(x, w),
    ge0_mld = log(mu) - wmean(log(x), w),
    ge1_theil = wmean((x / mu) * log(x / mu), w),
    ge2 = .5 * (wmean((x / mu)^2, w) - 1),
    atkinson_0_5 = 1 - ede05 / mu,
    atkinson_1 = 1 - ede1 / mu,
    atkinson_2 = 1 - ede2 / mu
  )
}

as_result <- function(name, x, w, baseline, poverty_line) {
  values <- indicators(x, w, poverty_line)
  data.frame(
    approach = name,
    indicator = names(values),
    value = unname(values),
    percent_change_vs_unadjusted = 100 * (unname(values) / unname(baseline) - 1),
    row.names = NULL
  )
}

# ---- Read and prepare data -------------------------------------------------
d <- read_dta(input_file)
valid <- is.finite(d$welfare) & d$welfare > 0 &
  is.finite(d$weight) & d$weight > 0 &
  is.finite(d$cpi2021) & d$cpi2021 > 0 &
  is.finite(d$icp2021) & d$icp2021 > 0

x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

# Same approximate $2.15 (2017 PPP) line expressed in 2021 PPP as prior work.
poverty_line <- 2.15 * (0.28 / 0.25)
baseline <- indicators(x, w, poverty_line)

# LIS upper fence on log welfare.
q1_log <- wq(log(x), w, .25)
q3_log <- wq(log(x), w, .75)
lis_ceiling <- exp(q3_log + 3 * (q3_log - q1_log))

# Classical 0.5% upper-tail adjustment, using a weighted cutoff.
top_share <- .005
top_cutoff <- wq(x, w, 1 - top_share)
top_records <- x > top_cutoff

# EU-SILC paper's upper-tail threshold rule:
# yh = max(min(2.5 * mean, Q(.98)), Q(.97)).
q97 <- wq(x, w, .97)
q98 <- wq(x, w, .98)
pareto_threshold <- max(min(2.5 * wmean(x, w), q98), q97)

# Robust Pareto fit (PDC) to the upper tail.  alpha=.005 defines observations
# that lie above the fitted 99.5th percentile of that tail as Pareto outliers.
pareto <- paretoTail(
  x,
  x0 = pareto_threshold,
  method = "thetaPDC",
  w = w,
  alpha = top_share
)
pareto_ceiling <- laeken:::qpareto(1 - top_share, pareto$x0, pareto$theta)
pareto_outliers <- rep(FALSE, length(x))
pareto_outliers[pareto$out] <- TRUE

# Model-based imputation: retain the 0.5% highest records but replace them
# using draws conditional on being above the observed top-0.5% threshold.
# The reported estimate is the average over 100 imputations; this removes
# incidental variation from one random draw while retaining the paper's logic.
set.seed(20260818)
imputations <- 100
imputed_values <- replicate(imputations, {
  x_imp <- x
  n_top <- sum(top_records)
  draws <- top_cutoff / runif(n_top)^(1 / pareto$theta)
  x_imp[top_records] <- sort(draws)
  indicators(x_imp, w, poverty_line)
})
imputation_mean <- rowMeans(imputed_values)
imputation_sd <- apply(imputed_values, 1, sd)

# ---- Compare approaches ----------------------------------------------------
x_lis <- pmin(x, lis_ceiling)
x_winsor <- pmin(x, top_cutoff)
x_pareto_shrink <- x
x_pareto_shrink[pareto_outliers] <- pareto_ceiling

results <- rbind(
  as_result("Unadjusted survey", x, w, baseline, poverty_line),
  as_result("Trim highest weighted 0.5%", x[!top_records], w[!top_records], baseline, poverty_line),
  as_result("Winsorize highest weighted 0.5%", x_winsor, w, baseline, poverty_line),
  as_result("LIS top-code: log-IQR k = 3", x_lis, w, baseline, poverty_line),
  as_result("Pareto PDC: shrink model-flagged values", x_pareto_shrink, w, baseline, poverty_line),
  data.frame(
    approach = "Pareto PDC: multiple imputation of top 0.5%",
    indicator = names(imputation_mean),
    value = unname(imputation_mean),
    percent_change_vs_unadjusted = 100 * (unname(imputation_mean) / unname(baseline) - 1),
    row.names = NULL
  )
)

imputation_uncertainty <- data.frame(
  approach = "Pareto PDC: multiple imputation of top 0.5%",
  indicator = names(imputation_sd),
  simulation_sd = unname(imputation_sd),
  imputations = imputations,
  row.names = NULL
)

method_notes <- data.frame(
  approach = c(
    "Unadjusted survey",
    "Trim highest weighted 0.5%",
    "Winsorize highest weighted 0.5%",
    "LIS top-code: log-IQR k = 3",
    "Pareto PDC: shrink model-flagged values",
    "Pareto PDC: multiple imputation of top 0.5%"
  ),
  treatment = c(
    "No adjustment.",
    "Removes records above the weighted 99.5th percentile.",
    "Caps records above the weighted 99.5th percentile at that percentile.",
    "Caps every value above exp(Q3(log y) + 3*IQR(log y)) at the LIS ceiling.",
    "Fits the upper tail with robust PDC Pareto; caps only values above its fitted 99.5th-tail percentile.",
    "Fits the same Pareto; replaces the highest weighted 0.5% with plausible Pareto draws, averaged over 100 replicates."
  ),
  interpretation = c(
    "Benchmark.",
    "Strong diagnostic: affected records receive zero influence.",
    "Classical top-coding with a fixed population share.",
    "LIS rule: transparent distribution-based ceiling.",
    "Model-based cap for values incompatible with the fitted tail.",
    "Model-based alternative that preserves a variable, plausible upper tail."
  ),
  stringsAsFactors = FALSE
)

settings <- data.frame(
  survey = "MWI_1997",
  valid_records = length(x),
  total_survey_weight = sum(w),
  poverty_line_ppp_day = poverty_line,
  weighted_p995_cutoff_ppp_day = top_cutoff,
  lis_ceiling_ppp_day = lis_ceiling,
  pareto_threshold_ppp_day = pareto_threshold,
  pareto_tail_records = pareto$k,
  pareto_pdc_alpha = pareto$theta,
  pareto_995_ceiling_ppp_day = pareto_ceiling,
  pareto_flagged_records = sum(pareto_outliers),
  pareto_flagged_weighted_share_percent = 100 * sum(w[pareto_outliers]) / sum(w),
  stringsAsFactors = FALSE
)

write.csv(results, file.path(output_dir, "mwi_1997_eusilc_approaches_comparison.csv"), row.names = FALSE)
write.csv(imputation_uncertainty, file.path(output_dir, "mwi_1997_eusilc_imputation_uncertainty.csv"), row.names = FALSE)
write.csv(method_notes, file.path(output_dir, "mwi_1997_eusilc_method_notes.csv"), row.names = FALSE)
write.csv(settings, file.path(output_dir, "mwi_1997_eusilc_settings.csv"), row.names = FALSE)

print(settings, row.names = FALSE)
print(results, row.names = FALSE)
message("Wrote EU-SILC comparison outputs to ", output_dir)
