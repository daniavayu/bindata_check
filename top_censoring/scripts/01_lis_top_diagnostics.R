#!/usr/bin/env Rscript

# Top-tail diagnostics for observed survey microdata.
# The script never changes input files. It writes summary CSVs only.

suppressPackageStartupMessages(library(haven))

input_dir <- file.path("01-input", "country")
target_file <- "MWI_1997.dta"
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

weighted_mean <- function(x, w) sum(x * w) / sum(w)

# Gini from the weighted Lorenz curve. Valid for non-negative welfare values.
weighted_gini <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & x >= 0 & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x) || sum(x * w) <= 0) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
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
  if (all(c("welfare", "cpi2021", "icp2021") %in% names(d))) {
    denominator <- d$cpi2021 * d$icp2021 * 365
    if (all(is.finite(denominator)) && all(denominator > 0)) {
      return(d$welfare / denominator)
    }
  }
  d$welfare
}

files <- file.path(input_dir, target_file)
if (!file.exists(files)) stop("Target survey file not found: ", files)

inventory <- list()
diagnostics <- list()

for (path in files) {
  d <- read_dta(path)
  file_name <- basename(path)
  columns <- names(d)
  grouped <- "bins" %in% columns
  weight_col <- infer_weight(columns)
  welfare_available <- "welfare" %in% columns

  inventory[[length(inventory) + 1]] <- data.frame(
    file = file_name,
    input_type = if (grouped) "grouped_distribution" else "microdata",
    rows = nrow(d),
    welfare_column = if (welfare_available) "welfare" else NA_character_,
    weight_column = weight_col,
    survey_year = if ("year" %in% columns) as.character(d$year[1]) else NA_character_,
    survey_name = if ("survname" %in% columns) as.character(d$survname[1]) else NA_character_,
    stringsAsFactors = FALSE
  )

  if (!welfare_available || is.na(weight_col)) next

  welfare <- to_ppp_day(d)
  weight <- as.numeric(d[[weight_col]])
  valid <- is.finite(welfare) & welfare > 0 & is.finite(weight) & weight > 0
  welfare <- welfare[valid]
  weight <- weight[valid]
  if (!length(welfare)) next

  log_quartiles <- weighted_quantile(log(welfare), weight, c(0.25, 0.75))
  log_iqr <- log_quartiles[2] - log_quartiles[1]

  for (k in c(2, 3, 4)) {
    ceiling <- exp(log_quartiles[2] + k * log_iqr)
    capped <- pmin(welfare, ceiling)
    top <- welfare > ceiling
    total_welfare <- sum(welfare * weight)

    diagnostics[[length(diagnostics) + 1]] <- data.frame(
      file = file_name,
      input_type = if (grouped) "grouped_distribution" else "microdata",
      survey_year = if ("year" %in% columns) as.character(d$year[1]) else NA_character_,
      survey_name = if ("survname" %in% columns) as.character(d$survname[1]) else NA_character_,
      weight_column = weight_col,
      lis_k = k,
      observations_positive_weight = length(welfare),
      population_weight = sum(weight),
      q1_log_welfare = log_quartiles[1],
      q3_log_welfare = log_quartiles[2],
      lis_ceiling_ppp_day = ceiling,
      weighted_share_above_ceiling = sum(weight[top]) / sum(weight),
      welfare_share_above_ceiling = sum(welfare[top] * weight[top]) / total_welfare,
      maximum_to_ceiling = max(welfare) / ceiling,
      mean_uncensored = weighted_mean(welfare, weight),
      mean_top_coded = weighted_mean(capped, weight),
      mean_change_pct = 100 * (weighted_mean(capped, weight) / weighted_mean(welfare, weight) - 1),
      gini_uncensored = weighted_gini(welfare, weight),
      gini_top_coded = weighted_gini(capped, weight),
      gini_change = weighted_gini(capped, weight) - weighted_gini(welfare, weight),
      stringsAsFactors = FALSE
    )
  }
}

write.csv(do.call(rbind, inventory), file.path(output_dir, "survey_inventory.csv"), row.names = FALSE)
if (length(diagnostics)) {
  write.csv(do.call(rbind, diagnostics), file.path(output_dir, "lis_top_diagnostics.csv"), row.names = FALSE)
}

message("Wrote outputs to ", output_dir)
