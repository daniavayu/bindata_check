#!/usr/bin/env Rscript

# Pareto tail sensitivity for Malawi 1997: original data versus the same p90
# tail after excluding the two largest welfare observations.

suppressPackageStartupMessages(library(haven))

d <- read_dta(file.path("01-input", "country", "MWI_1997.dta"))
output_dir <- file.path("top_censoring", "outputs", "pareto")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

valid <- is.finite(d$welfare) & d$welfare > 0 & is.finite(d$weight) & d$weight > 0
x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

# Keep the original p90 threshold fixed in both fits for comparability.
threshold <- wq(x, w, .90)
top_two <- order(x, decreasing = TRUE)[1:2]
keep <- rep(TRUE, length(x)); keep[top_two] <- FALSE

fit_pareto <- function(x, w, threshold) {
  tail <- x >= threshold
  tx <- x[tail]; tw <- w[tail]
  alpha <- sum(tw) / sum(tw * log(tx / threshold))
  o <- order(tx); tx <- tx[o]; tw <- tw[o]
  list(
    alpha = alpha,
    records = length(tx),
    x = tx,
    empirical = rev(cumsum(rev(tw))) / sum(tw),
    fitted = (threshold / tx) ^ alpha
  )
}

all_fit <- fit_pareto(x, w, threshold)
trimmed_fit <- fit_pareto(x[keep], w[keep], threshold)

summary <- data.frame(
  scenario = c("all_records", "excluding_two_largest"),
  p90_threshold_ppp_day = threshold,
  tail_records = c(all_fit$records, trimmed_fit$records),
  pareto_alpha = c(all_fit$alpha, trimmed_fit$alpha)
)
write.csv(summary, file.path(output_dir, "mwi_1997_pareto_without_two_summary.csv"), row.names = FALSE)

png(file.path(output_dir, "mwi_1997_pareto_without_two.png"), width = 2200, height = 1100, res = 180)
par(mfrow = c(1, 2))
for (item in list(list(fit = all_fit, title = "All records"), list(fit = trimmed_fit, title = "Excluding two largest values"))) {
  plot(item$fit$x, item$fit$empirical, log = "xy", type = "s", lwd = 1.5, col = "#0072B2",
       xlab = "Welfare (2021 PPP USD/day)", ylab = "Share of p90 tail at or above welfare",
       main = paste("Malawi 1997:", item$title))
  lines(item$fit$x, item$fit$fitted, lwd = 1.5, col = "#D55E00")
  legend("topright", c("Observed tail", "Pareto"), col = c("#0072B2", "#D55E00"), lty = 1, lwd = 1.5, bty = "n")
}
dev.off()

print(summary, row.names = FALSE)
