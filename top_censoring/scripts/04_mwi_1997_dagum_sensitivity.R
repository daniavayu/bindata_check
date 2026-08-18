#!/usr/bin/env Rscript

# Dagum Type I fits for Malawi 1997, with and without the two largest welfare
# values. This is a distributional sensitivity check, not a data-cleaning rule.

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

# Dagum Type I: F(x) = [1 + (x / b)^(-a)]^(-p), with a, b, p > 0.
dagum_cdf <- function(x, a, b, p) {
  z <- -a * (log(x) - log(b))
  exp(-p * log1pexp(z))
}

fit_dagum <- function(x, w) {
  normalized_w <- w / sum(w)
  log_x <- log(x)
  log_b_start <- log(weighted_quantile(x, w, 0.50))

  negative_log_likelihood <- function(theta) {
    a <- exp(theta[1])
    b <- exp(theta[2])
    p <- exp(theta[3])
    z <- -a * (log_x - log(b))
    log_density <- log(a) + log(p) - log_x + z - (p + 1) * log1pexp(z)
    -sum(normalized_w * log_density)
  }

  starts <- expand.grid(a = c(0.5, 1, 2), p = c(0.5, 1, 2))
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = c(log(starts$a[i]), log_b_start, log(starts$p[i])),
      fn = negative_log_likelihood,
      method = "L-BFGS-B",
      lower = c(log(0.01), log(min(x)) - 2, log(0.01)),
      upper = c(log(100), log(max(x)) + 2, log(100))
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  params <- exp(best$par)
  names(params) <- c("a", "b", "p")

  order_index <- order(x)
  empirical_cdf <- cumsum(w[order_index]) / sum(w)
  fitted_cdf <- dagum_cdf(x[order_index], params["a"], params["b"], params["p"])

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

fit_all <- fit_dagum(x, w)
fit_trimmed <- fit_dagum(x[without_top_two], w[without_top_two])

summary <- data.frame(
  scenario = c("all_records", "excluding_two_largest"),
  records = c(length(x), sum(without_top_two)),
  population_weight = c(sum(w), sum(w[without_top_two])),
  dagum_a = c(fit_all$parameters["a"], fit_trimmed$parameters["a"]),
  dagum_b_ppp_day = c(fit_all$parameters["b"], fit_trimmed$parameters["b"]),
  dagum_p = c(fit_all$parameters["p"], fit_trimmed$parameters["p"]),
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
  file.path(output_dir, "mwi_1997_dagum_sensitivity_summary.csv"),
  row.names = FALSE
)

plot_fit <- function(x, w, fit, title) {
  order_index <- order(x)
  sorted_x <- x[order_index]
  sorted_w <- w[order_index]
  empirical_survival <- 1 - cumsum(sorted_w) / sum(sorted_w) + sorted_w / sum(sorted_w)
  fitted_survival <- 1 - dagum_cdf(
    sorted_x,
    fit$parameters["a"],
    fit$parameters["b"],
    fit$parameters["p"]
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
    col = "#D55E00"
  )
  legend(
    "topright",
    legend = c("Observed weighted distribution", "Fitted Dagum Type I"),
    col = c("#0072B2", "#D55E00"),
    lty = 1,
    lwd = 1.5,
    bty = "n"
  )
}

png(
  file.path(output_dir, "mwi_1997_dagum_sensitivity.png"),
  width = 2200,
  height = 1100,
  res = 180
)
par(mfrow = c(1, 2))
plot_fit(x, w, fit_all, "Malawi 1997: Dagum fit, all records")
plot_fit(
  x[without_top_two],
  w[without_top_two],
  fit_trimmed,
  "Malawi 1997: Dagum fit, excluding two largest values"
)
dev.off()

print(summary, row.names = FALSE)
message("Wrote Dagum outputs to ", output_dir)
