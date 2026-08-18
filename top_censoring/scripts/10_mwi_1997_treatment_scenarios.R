#!/usr/bin/env Rscript

# Compare alternative treatments of Malawi 1997 upper-tail observations.

suppressPackageStartupMessages(library(haven))

input_file <- file.path("01-input", "country", "MWI_1997.dta")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_dta(input_file)
valid <- is.finite(d$welfare) & d$welfare > 0 &
  is.finite(d$weight) & d$weight > 0
x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}
wmean <- function(x, w) sum(x * w) / sum(w)
wgini <- function(x, w) {
  o <- order(x); x <- x[o]; w <- w[o]
  pop <- c(0, cumsum(w) / sum(w)); inc <- c(0, cumsum(x * w) / sum(x * w))
  1 - sum(diff(pop) * (head(inc, -1) + tail(inc, -1)))
}
indicators <- function(x, w) {
  mu <- wmean(x, w)
  ede2 <- 1 / wmean(1 / x, w)
  c(mean_ppp_day = mu, gini = wgini(x, w),
    mld = log(mu) - wmean(log(x), w), atkinson_2 = 1 - ede2 / mu)
}

q1 <- wq(log(x), w, .25); q3 <- wq(log(x), w, .75)
lis_ceiling <- exp(q3 + 3 * (q3 - q1))
top_two <- order(x, decreasing = TRUE)[1:2]
above_lis <- x > lis_ceiling

x_cap_two <- x
x_cap_two[top_two] <- lis_ceiling
x_cap_all_lis <- pmin(x, lis_ceiling)
keep_without_two <- rep(TRUE, length(x)); keep_without_two[top_two] <- FALSE

scenarios <- list(
  original = list(x = x, w = w),
  cap_two_at_lis_ceiling = list(x = x_cap_two, w = w),
  cap_all_lis_values = list(x = x_cap_all_lis, w = w),
  delete_two_diagnostic_only = list(x = x[keep_without_two], w = w[keep_without_two])
)

baseline <- indicators(x, w)
result <- do.call(rbind, lapply(names(scenarios), function(name) {
  values <- indicators(scenarios[[name]]$x, scenarios[[name]]$w)
  data.frame(
    scenario = name,
    indicator = names(values),
    value = unname(values),
    percent_change_vs_original = 100 * (unname(values) / unname(baseline) - 1)
  )
}))

write.csv(result, file.path(output_dir, "mwi_1997_treatment_scenarios.csv"), row.names = FALSE)
print(result, row.names = FALSE)
message("LIS ceiling: ", round(lis_ceiling, 2), " PPP 2021 USD/day")
