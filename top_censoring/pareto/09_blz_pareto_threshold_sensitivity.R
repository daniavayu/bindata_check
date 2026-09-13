#!/usr/bin/env Rscript

# Pareto threshold sensitivity for Belize 1993 and 1997.

suppressPackageStartupMessages(library(haven))

input_dir <- file.path("01-input", "country")
output_dir <- file.path("top_censoring", "pareto", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

fit_pareto <- function(x, w, percentile) {
  threshold <- wq(x, w, percentile)
  in_tail <- x >= threshold
  tx <- x[in_tail]; tw <- w[in_tail]
  alpha <- sum(tw) / sum(tw * log(tx / threshold))
  o <- order(tx)
  tx <- tx[o]; tw <- tw[o]
  empirical <- rev(cumsum(rev(tw))) / sum(tw)
  fitted <- (threshold / tx) ^ alpha
  list(threshold = threshold, alpha = alpha, x = tx, empirical = empirical, fitted = fitted)
}

all_results <- list()
plots <- list()
for (file in c("BLZ_1993.dta", "BLZ_1997.dta")) {
  d <- read_dta(file.path(input_dir, file))
  valid <- is.finite(d$welfare) & d$welfare > 0 & is.finite(d$weight) & d$weight > 0
  x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
  w <- d$weight[valid]
  year <- as.character(d$year[1])
  q1 <- wq(log(x), w, .25); q3 <- wq(log(x), w, .75)
  lis <- exp(q3 + 3 * (q3 - q1))

  for (p in c(.90, .95, .99)) {
    fit <- fit_pareto(x, w, p)
    above_lis <- fit$x > lis
    # The largest empirical/model survival ratio above the LIS ceiling.
    max_ratio <- if (any(above_lis)) max(fit$empirical[above_lis] / fit$fitted[above_lis]) else NA_real_
    all_results[[length(all_results) + 1]] <- data.frame(
      year = as.integer(year), tail_start = paste0("p", p * 100),
      threshold_ppp_day = fit$threshold, pareto_alpha = fit$alpha,
      maximum_empirical_to_pareto_above_lis = max_ratio
    )
    plots[[paste(year, p)]] <- fit
  }
}

summary <- do.call(rbind, all_results)
write.csv(summary, file.path(output_dir, "blz_pareto_threshold_sensitivity.csv"), row.names = FALSE)

png(file.path(output_dir, "blz_pareto_threshold_sensitivity.png"), width = 2400, height = 1500, res = 180)
par(mfrow = c(2, 3))
for (year in c("1993", "1997")) {
  for (p in c(.90, .95, .99)) {
    fit <- plots[[paste(year, p)]]
    plot(fit$x, fit$empirical, log = "xy", type = "s", lwd = 1.2, col = "#0072B2",
         xlab = "Welfare (2021 PPP USD/day)", ylab = "Tail survival share",
         main = paste("Belize", year, "- Pareto from p", p * 100, sep = ""))
    lines(fit$x, fit$fitted, col = "#D55E00", lwd = 1.2)
    legend("topright", c("Observed", "Pareto"), col = c("#0072B2", "#D55E00"),
           lty = 1, lwd = 1.2, bty = "n", cex = .8)
  }
}
dev.off()

print(summary, row.names = FALSE)
