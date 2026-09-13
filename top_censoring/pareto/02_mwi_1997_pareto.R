#!/usr/bin/env Rscript

# Malawi 1997 upper-tail Pareto diagnostic.
# Inputs are read only; figures and summaries are written to top_censoring/pareto/outputs.

suppressPackageStartupMessages(library(haven))

input_file <- file.path("01-input", "country", "MWI_1997.dta")
output_dir <- file.path("top_censoring", "pareto", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_dta(input_file)

valid <- is.finite(d$welfare) & d$welfare > 0 &
  is.finite(d$weight) & d$weight > 0

# Convert to 2021 PPP USD per person per day.
x <- d$welfare[valid] /
  (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

weighted_quantile <- function(x, w, p) {
  order_index <- order(x)
  sorted_x <- x[order_index]
  sorted_w <- w[order_index]
  sorted_x[which(cumsum(sorted_w) / sum(sorted_w) >= p)[1]]
}

# Fit a Pareto distribution to the upper 10 percent of the weighted distribution.
tail_percentile <- 0.90
threshold <- weighted_quantile(x, w, tail_percentile)
in_tail <- x >= threshold
tail_x <- x[in_tail]
tail_w <- w[in_tail]

# Weighted maximum-likelihood estimator of the Pareto shape parameter.
pareto_alpha <- sum(tail_w) / sum(tail_w * log(tail_x / threshold))

# LIS ceiling: three IQRs above Q3 of log welfare.
q1_log <- weighted_quantile(log(x), w, 0.25)
q3_log <- weighted_quantile(log(x), w, 0.75)
lis_ceiling <- exp(q3_log + 3 * (q3_log - q1_log))

# Weighted empirical survival curve and fitted Pareto survival curve.
order_index <- order(tail_x)
tail_x <- tail_x[order_index]
tail_w <- tail_w[order_index]
empirical_survival <- rev(cumsum(rev(tail_w))) / sum(tail_w)
pareto_survival <- (threshold / tail_x)^pareto_alpha

png(
  file.path(output_dir, "mwi_1997_pareto_p90_ccdf.png"),
  width = 1800,
  height = 1200,
  res = 180
)

plot(
  tail_x,
  empirical_survival,
  log = "xy",
  type = "s",
  lwd = 2,
  col = "#0072B2",
  xlab = "Welfare (2021 PPP USD per person per day, log scale)",
  ylab = "Share of p90 tail at or above welfare (log scale)",
  main = "Malawi 1997: observed upper tail versus Pareto fit"
)
lines(tail_x, pareto_survival, lwd = 2, col = "#D55E00")
abline(v = lis_ceiling, col = "#009E73", lty = 2, lwd = 2)
legend(
  "topright",
  legend = c(
    "Observed weighted tail",
    "Pareto fitted from p90",
    "LIS k = 3 ceiling"
  ),
  col = c("#0072B2", "#D55E00", "#009E73"),
  lty = c(1, 1, 2),
  lwd = 2,
  bty = "n"
)
dev.off()

summary <- data.frame(
  survey = "MWI_1997",
  tail_percentile = tail_percentile,
  threshold_ppp_day = threshold,
  pareto_alpha = pareto_alpha,
  lis_k = 3,
  lis_ceiling_ppp_day = lis_ceiling,
  stringsAsFactors = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "mwi_1997_pareto_p90_summary.csv"),
  row.names = FALSE
)

message("Wrote Pareto outputs to ", output_dir)
