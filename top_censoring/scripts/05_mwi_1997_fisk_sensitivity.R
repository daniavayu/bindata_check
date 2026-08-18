#!/usr/bin/env Rscript

# Fisk (log-logistic) fits for Malawi 1997, with and without the two largest
# welfare values. Fisk is a simpler full-distribution benchmark than Dagum.

suppressPackageStartupMessages(library(haven))

input_file <- file.path("01-input", "country", "MWI_1997.dta")
output_dir <- file.path("top_censoring", "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_dta(input_file)
valid <- is.finite(d$welfare) & d$welfare > 0 &
  is.finite(d$weight) & d$weight > 0
x <- d$welfare[valid] / (d$cpi2021[valid] * d$icp2021[valid] * 365)
w <- d$weight[valid]

weighted_quantile <- function(x, w, p) {
  order_index <- order(x)
  sorted_x <- x[order_index]
  sorted_w <- w[order_index]
  sorted_x[which(cumsum(sorted_w) / sum(sorted_w) >= p)[1]]
}

log1pexp <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))

# Fisk CDF: F(x) = 1 / [1 + (x / b)^(-a)], with a, b > 0.
fisk_cdf <- function(x, a, b) {
  z <- -a * (log(x) - log(b))
  exp(-log1pexp(z))
}

fit_fisk <- function(x, w) {
  normalized_w <- w / sum(w)
  log_x <- log(x)
  log_b_start <- log(weighted_quantile(x, w, 0.50))

  negative_log_likelihood <- function(theta) {
    a <- exp(theta[1])
    b <- exp(theta[2])
    z <- -a * (log_x - log(b))
    log_density <- log(a) - log_x + z - 2 * log1pexp(z)
    -sum(normalized_w * log_density)
  }

  starts <- c(0.5, 1, 2, 4)
  fits <- lapply(starts, function(a_start) {
    optim(
      par = c(log(a_start), log_b_start),
      fn = negative_log_likelihood,
      method = "L-BFGS-B",
      lower = c(log(0.01), log(min(x)) - 2),
      upper = c(log(100), log(max(x)) + 2)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  params <- exp(best$par)
  names(params) <- c("a", "b")

  order_index <- order(x)
  empirical_cdf <- cumsum(w[order_index]) / sum(w)
  fitted_cdf <- fisk_cdf(x[order_index], params["a"], params["b"])

  list(
    parameters = params,
    weighted_ks = max(abs(empirical_cdf - fitted_cdf)),
    negative_log_likelihood = best$value,
    convergence = best$convergence
  )
}

top_two <- order(x, decreasing = TRUE)[1:2]
without_top_two <- rep(TRUE, length(x))
without_top_two[top_two] <- FALSE

fit_all <- fit_fisk(x, w)
fit_trimmed <- fit_fisk(x[without_top_two], w[without_top_two])

summary <- data.frame(
  scenario = c("all_records", "excluding_two_largest"),
  records = c(length(x), sum(without_top_two)),
  fisk_a = c(fit_all$parameters["a"], fit_trimmed$parameters["a"]),
  fisk_b_ppp_day = c(fit_all$parameters["b"], fit_trimmed$parameters["b"]),
  weighted_ks = c(fit_all$weighted_ks, fit_trimmed$weighted_ks),
  negative_log_likelihood = c(
    fit_all$negative_log_likelihood,
    fit_trimmed$negative_log_likelihood
  ),
  convergence_code = c(fit_all$convergence, fit_trimmed$convergence),
  row.names = NULL
)
write.csv(
  summary,
  file.path(output_dir, "mwi_1997_fisk_sensitivity_summary.csv"),
  row.names = FALSE
)

plot_fit <- function(x, w, fit, title) {
  order_index <- order(x)
  sorted_x <- x[order_index]
  sorted_w <- w[order_index]
  empirical_survival <- 1 - cumsum(sorted_w) / sum(sorted_w) + sorted_w / sum(sorted_w)
  fitted_survival <- 1 - fisk_cdf(
    sorted_x,
    fit$parameters["a"],
    fit$parameters["b"]
  )

  plot(
    sorted_x,
    pmax(empirical_survival, .Machine$double.xmin),
    log = "xy",
    type = "s",
    lwd = 1.5,
    col = "#0072B2",
    xlab = "Welfare (2021 PPP USD per person per day, log scale)",
    ylab = "Share at or above welfare (log scale)",
    main = title
  )
  lines(
    sorted_x,
    pmax(fitted_survival, .Machine$double.xmin),
    lwd = 1.5,
    col = "#CC79A7"
  )
  legend(
    "topright",
    legend = c("Observed weighted distribution", "Fitted Fisk"),
    col = c("#0072B2", "#CC79A7"),
    lty = 1,
    lwd = 1.5,
    bty = "n"
  )
}

png(
  file.path(output_dir, "mwi_1997_fisk_sensitivity.png"),
  width = 2200,
  height = 1100,
  res = 180
)
par(mfrow = c(1, 2))
plot_fit(x, w, fit_all, "Malawi 1997: Fisk fit, all records")
plot_fit(
  x[without_top_two],
  w[without_top_two],
  fit_trimmed,
  "Malawi 1997: Fisk fit, excluding two largest values"
)
dev.off()

print(summary, row.names = FALSE)
message("Wrote Fisk outputs to ", output_dir)
