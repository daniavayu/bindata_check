#!/usr/bin/env Rscript

# Sensitivity of Malawi 1997 welfare, poverty, and inequality indicators to
# removing the two largest welfare observations. This is a diagnostic only,
# not a recommendation to remove observations from the survey.

suppressPackageStartupMessages(library(haven))

input_file <- file.path("01-input", "country", "MWI_1997.dta")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_dta(input_file)
valid <- is.finite(d$welfare) & d$welfare > 0 &
  is.finite(d$weight) & d$weight > 0
x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

weighted_mean <- function(x, w) sum(x * w) / sum(w)

weighted_gini <- function(x, w) {
  order_index <- order(x)
  x <- x[order_index]
  w <- w[order_index]
  population <- c(0, cumsum(w) / sum(w))
  welfare <- c(0, cumsum(x * w) / sum(x * w))
  1 - sum(diff(population) * (head(welfare, -1) + tail(welfare, -1)))
}

indicators <- function(x, w, poverty_line) {
  mean_welfare <- weighted_mean(x, w)
  poor <- x < poverty_line
  gap <- numeric(length(x))
  gap[poor] <- (poverty_line - x[poor]) / poverty_line

  # The MLD and Atkinson(2) use the positive welfare support in this survey.
  mean_log_deviation <- log(mean_welfare) - weighted_mean(log(x), w)
  equally_distributed_equivalent <- 1 / weighted_mean(1 / x, w)
  atkinson_2 <- 1 - equally_distributed_equivalent / mean_welfare

  c(
    mean_ppp_day = mean_welfare,
    headcount_ratio = sum(w[poor]) / sum(w),
    fgt1 = weighted_mean(gap, w),
    fgt2 = weighted_mean(gap ^ 2, w),
    watts = weighted_mean(ifelse(poor, log(poverty_line / x), 0), w),
    gini = weighted_gini(x, w),
    mld = mean_log_deviation,
    atkinson_2 = atkinson_2
  )
}

# The paper's $2.15/day 2017-PPP line is approximately $2.41 in 2021 PPP,
# using the paper's stated $0.25 (2017) to $0.28 (2021) conversion.
poverty_line <- 2.15 * (0.28 / 0.25)
top_two <- order(x, decreasing = TRUE)[1:2]
keep_without_top_two <- rep(TRUE, length(x))
keep_without_top_two[top_two] <- FALSE

all_records <- indicators(x, w, poverty_line)
without_top_two <- indicators(
  x[keep_without_top_two],
  w[keep_without_top_two],
  poverty_line
)

result <- data.frame(
  indicator = names(all_records),
  all_records = unname(all_records),
  excluding_two_largest = unname(without_top_two),
  absolute_change = unname(without_top_two - all_records),
  percent_change = unname(100 * (without_top_two / all_records - 1)),
  row.names = NULL
)

write.csv(
  result,
  file.path(output_dir, "mwi_1997_two_extremes_index_sensitivity.csv"),
  row.names = FALSE
)

print(result, row.names = FALSE)
message("Poverty line: ", round(poverty_line, 2), " PPP 2021 USD/day")
