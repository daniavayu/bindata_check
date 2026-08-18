#!/usr/bin/env Rscript

# Initial top-tail analysis for Belize 1993 and Belize 1997.

suppressPackageStartupMessages(library(haven))

input_dir <- file.path("01-input", "country")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- c("BLZ_1993.dta", "BLZ_1997.dta")

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
indices <- function(x, w) {
  mu <- wmean(x, w)
  ede2 <- 1 / wmean(1 / x, w)
  c(mean_ppp_day = mu, gini = wgini(x, w),
    mld = log(mu) - wmean(log(x), w), atkinson_2 = 1 - ede2 / mu)
}

tail_results <- list()
index_results <- list()
plot_data <- list()

for (file in files) {
  d <- read_dta(file.path(input_dir, file))
  valid <- is.finite(d$welfare) & d$welfare > 0 & is.finite(d$weight) & d$weight > 0
  x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
  w <- d$weight[valid]
  year <- as.integer(d$year[1])

  q1 <- wq(log(x), w, .25); q3 <- wq(log(x), w, .75)
  lis <- exp(q3 + 3 * (q3 - q1))
  above <- x > lis

  u <- wq(x, w, .90)
  tail <- x >= u
  tx <- x[tail]; tw <- w[tail]
  alpha <- sum(tw) / sum(tw * log(tx / u))
  top_two <- order(x, decreasing = TRUE)[1:2]
  top_x <- x[top_two]
  empirical <- vapply(top_x, function(z) sum(tw[tx >= z]) / sum(tw), numeric(1))
  predicted <- (u / top_x) ^ alpha

  tail_results[[length(tail_results) + 1]] <- data.frame(
    country = "BLZ", year = year, records = length(x), lis_ceiling_ppp_day = lis,
    records_above_lis = sum(above), population_share_above_lis_pct = 100 * sum(w[above]) / sum(w),
    welfare_share_above_lis_pct = 100 * sum(x[above] * w[above]) / sum(x * w),
    pareto_p90_threshold = u, pareto_alpha = alpha,
    top1_ppp_day = top_x[1], top1_empirical_to_pareto = empirical[1] / predicted[1],
    top2_ppp_day = top_x[2], top2_empirical_to_pareto = empirical[2] / predicted[2]
  )

  keep <- rep(TRUE, length(x)); keep[top_two] <- FALSE
  all_values <- indices(x, w)
  trimmed_values <- indices(x[keep], w[keep])
  index_results[[length(index_results) + 1]] <- data.frame(
    country = "BLZ", year = year, indicator = names(all_values),
    all_records = unname(all_values), excluding_two_largest = unname(trimmed_values),
    percent_change = 100 * (unname(trimmed_values) / unname(all_values) - 1)
  )

  o <- order(tx)
  plot_data[[as.character(year)]] <- list(
    x = tx[o], empirical = rev(cumsum(rev(tw[o]))) / sum(tw),
    pareto = (u / tx[o]) ^ alpha, lis = lis
  )
}

tail_summary <- do.call(rbind, tail_results)
index_summary <- do.call(rbind, index_results)
write.csv(tail_summary, file.path(output_dir, "blz_1993_1997_pareto_summary.csv"), row.names = FALSE)
write.csv(index_summary, file.path(output_dir, "blz_1993_1997_index_sensitivity.csv"), row.names = FALSE)

png(file.path(output_dir, "blz_1993_1997_pareto_ccdf.png"), width = 2200, height = 1100, res = 180)
par(mfrow = c(1, 2))
for (year in c("1993", "1997")) {
  z <- plot_data[[year]]
  plot(z$x, z$empirical, log = "xy", type = "s", lwd = 1.5, col = "#0072B2",
       xlab = "Welfare (2021 PPP USD/day, log scale)",
       ylab = "Share of p90 tail at or above welfare (log scale)",
       main = paste("Belize", year, ": Pareto fit from p90"))
  lines(z$x, z$pareto, lwd = 1.5, col = "#D55E00")
  abline(v = z$lis, lty = 2, lwd = 1.5, col = "#009E73")
  legend("topright", c("Observed tail", "Pareto", "LIS ceiling"),
         col = c("#0072B2", "#D55E00", "#009E73"), lty = c(1, 1, 2), lwd = 1.5, bty = "n")
}
dev.off()

print(tail_summary, row.names = FALSE)
print(index_summary, row.names = FALSE)
