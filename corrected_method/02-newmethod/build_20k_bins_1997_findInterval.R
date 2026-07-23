library(haven)
library(data.table)
library(fastverse)
library(openxlsx)

fastverse_conflicts()
conflicted::conflicts_prefer(data.table::between())
conflicted::conflicts_prefer(collapse::funique())
conflicted::conflicts_prefer(data.table::first())

# =============================================================================
# Builds 20,000 bins for MWI reporting year 1997 using the generalized-
# Lorenz-curve approach (preserves the exact weighted mean, unlike PIP's
# quantile-point method which only approximates it).
#
# This is a drop-in replacement for build_20k_bins_1997.R: same inputs, same
# outputs, only the interpolation step (step 4) changes.
#
#   - Original: base R approx(..., method = "linear")
#   - Here:     manual vectorized linear interpolation via findInterval()
#
# Why: findInterval() gives ~2x the speed of approx() on this workload (see
# benchmark_lorenz_interpolation.R), which matters once this runs across
# every country/reporting-year in the pipeline, not just MWI 1997. Accuracy
# was verified in that benchmark: max abs difference vs. approx() ~ 1e-9.
# =============================================================================

# ---- Paths ------------------------------------------------------------------

ROOT   <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
INPUT  <- file.path(ROOT, "corrected_method", "01-processing")
OUTPUT <- file.path(ROOT, "corrected_method", "03-outputs")
mwi_path <- file.path(INPUT, "mwi_1997_reportingyear_processed.dta")

# ---- Helper: Lorenz-curve interpolation (the findInterval upgrade) ----------

#' Linearly interpolate cumulative welfare at a set of population-share
#' boundaries, using findInterval() instead of approx().
#'
#' @param p_cum    Strictly increasing cumulative population share, in [0, 1].
#' @param cum_inc  Cumulative (weighted) welfare at each point in p_cum.
#' @param p_bounds Population-share points at which to interpolate cum_inc.
#' @return Numeric vector, same length as p_bounds.
lorenz_interp <- function(p_cum, cum_inc, p_bounds) {
  n <- length(p_cum)
  idx <- findInterval(p_bounds, p_cum, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1L), n - 1L)

  x0 <- p_cum[idx];   x1 <- p_cum[idx + 1L]
  y0 <- cum_inc[idx]; y1 <- cum_inc[idx + 1L]

  y0 + (y1 - y0) * (p_bounds - x0) / (x1 - x0)
}

# ---- Main function ----------------------------------------------------------

#' Build nobs equal-population bins per reporting_level, preserving the exact
#' weighted mean via the generalized Lorenz curve.
#'
#' @param df   data.table with columns welfare, weight, reporting_level.
#' @param nobs Number of bins to build (default 20,000).
#' @return data.table with columns welfare, weight, reporting_level.
get_refy_quantiles_mean <- function(df, nobs = 2e4) {

  setorder(df, welfare)
  df_attr <- attributes(df)
  rls <- funique(df$reporting_level)

  # 1. Define the nobs + 1 exact boundaries of the bins (from 0 to 1)
  p_bounds <- seq(0, 1, length.out = nobs + 1)

  qx <- lapply(rls, \(rl) {
    # Filter out 0-weight rows to avoid breaking the cumulative interpolation
    x <- df[reporting_level == rl & weight > 0]
    xpop <- fsum(x$weight)

    # 2. Cumulative population share (p_cum)
    p_cum <- c(0, fcumsum(x$weight) / xpop)
    p_cum[length(p_cum)] <- 1  # guard against floating-point rounding (e.g. 0.999999...)

    # 3. Cumulative total welfare (Generalized Lorenz curve)
    cum_inc <- c(0, fcumsum(x$weight * x$welfare))

    # 4. Interpolate the exact cumulative welfare at the bin boundaries.
    #    Linear interpolation perfectly splits households that overlap two bins.
    interpolated_inc <- lorenz_interp(p_cum, cum_inc, p_bounds)

    # 5. Total welfare inside each bin is the difference between boundaries
    bin_inc <- diff(interpolated_inc)

    # 6. The mean welfare of the bin is its total welfare / its population size
    bin_welfare <- bin_inc / (xpop / nobs)

    data.table(welfare = bin_welfare, weight = xpop / nobs, reporting_level = rl)
  }) |>
    rowbind()

  # Restore attributes
  for (nm in names(df_attr)) {
    if (!nm %in% c("dim", "row.names", "names", "class", ".internal.selfref")) {
      setattr(qx, nm, df_attr[[nm]])
    }
  }

  qx
}

# ---- Run ----------------------------------------------------------------

my_data <- as.data.table(read_dta(mwi_path))
my_data[, reporting_level := "national"]

cat("Input data (reporting year 1997):\n")
cat("  Rows:", nrow(my_data), "\n")
cat("  Weighted mean:", weighted.mean(my_data$welfare, my_data$weight), "\n\n")

cat("Creating 20k bins (findInterval upgrade)...\n")
result <- get_refy_quantiles_mean(my_data, nobs = 2e4)

cat("\nResult:\n")
cat("  Rows:", nrow(result), "\n")
cat("  Weighted mean:", weighted.mean(result$welfare, result$weight), "\n")
cat("  Total weight:", sum(result$weight), "\n")

# ---- Save -----------------------------------------------------------------

output_path <- file.path(OUTPUT, "MWI_1997_reportingyear_20k_mean_findInterval.xlsx")
write.xlsx(result, output_path)
cat("\nSaved to:", output_path, "\n")
