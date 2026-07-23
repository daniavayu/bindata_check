library(haven)
library(data.table)
library(fastverse)

fastverse_conflicts()
conflicted::conflicts_prefer(data.table::between())
conflicted::conflicts_prefer(collapse::funique())
conflicted::conflicts_prefer(data.table::first())

# Benchmark: which implementation of the Lorenz-curve interpolation step
# (the core of get_refy_quantiles_mean(), see build_20k_bins_1997.R) is
# fastest? Accuracy was already validated (weighted mean matches exactly);
# this only compares speed, since the same interpolation has to run for
# every country/reporting-year in the pipeline.

ROOT <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
mwi_path <- file.path(ROOT, "corrected_method", "01-processing", "mwi_1997_reportingyear_processed.dta")

my_data <- as.data.table(read_dta(mwi_path))
setorder(my_data, welfare)
x <- my_data[weight > 0]

nobs <- 2e4
p_bounds <- seq(0, 1, length.out = nobs + 1)
xpop <- fsum(x$weight)

p_cum <- c(0, fcumsum(x$weight) / xpop)
p_cum[length(p_cum)] <- 1
cum_inc <- c(0, fcumsum(x$weight * x$welfare))

# --- Three candidate implementations of the interpolation step ---------

# 1. Current implementation: base R approx()
interp_approx <- function() {
  approx(x = p_cum, y = cum_inc, xout = p_bounds, method = "linear")$y
}

# 2. Manual vectorized linear interpolation using findInterval()
interp_findInterval <- function() {
  idx <- findInterval(p_bounds, p_cum, rightmost.closed = TRUE, all.inside = TRUE)
  denom <- p_cum[idx + 1L] - p_cum[idx]
  frac <- fifelse(denom == 0, 0, (p_bounds - p_cum[idx]) / denom)
  cum_inc[idx] + frac * (cum_inc[idx + 1L] - cum_inc[idx])
}

# 3. stats::splinefun with method = "linear" (builds an interpolating
#    function once, then evaluates it)
interp_splinefun <- function() {
  splinefun(x = p_cum, y = cum_inc, method = "linear")(p_bounds)
}

# --- Sanity check: all three must give (numerically) the same bins -----
r1 <- interp_approx()
r2 <- interp_findInterval()
r3 <- interp_splinefun()
stopifnot(
  isTRUE(all.equal(r1, r2)),
  isTRUE(all.equal(r1, r3))
)
cat("Sanity check passed: all three approaches yield the same result.\n\n")

# --- Benchmark -----------------------------------------------------------
bench <- microbenchmark::microbenchmark(
  times = 100,
  approx = {
    interp_approx()
  },
  findInterval = {
    interp_findInterval()
  },
  splinefun = {
    interp_splinefun()
  }
)

print(bench)

if (requireNamespace("highcharter", quietly = TRUE)) {
  hc_dt <- highcharter::data_to_boxplot(bench,
                                        time,
                                        expr,
                                        add_outliers = FALSE,
                                        name = "Time in milliseconds")

  highcharter::highchart() |>
    highcharter::hc_xAxis(type = "category") |>
    highcharter::hc_chart(inverted = TRUE) |>
    highcharter::hc_add_series_list(hc_dt)

} else {
  boxplot(bench, outline = FALSE)
}
