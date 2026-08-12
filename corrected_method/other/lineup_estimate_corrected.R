# Corrected version of get_refy_quantiles() from:
# PIP-Technical-Team/pip_ingestion_pipeline/R/lineup_estimate.R
#
# This function keeps the production function name, arguments, output
# columns, bin weights, ordering, and attributes. The welfare value in each
# bin is calculated from the generalized Lorenz curve so the weighted mean is
# preserved exactly.

lorenz_interp <- function(p_cum, cum_inc, p_bounds) {
  n <- length(p_cum)
  idx <- findInterval(p_bounds, p_cum, rightmost.closed = TRUE)
  idx <- pmin(pmax(idx, 1L), n - 1L)

  x0 <- p_cum[idx]
  x1 <- p_cum[idx + 1L]
  y0 <- cum_inc[idx]
  y1 <- cum_inc[idx + 1L]

  y0 + (y1 - y0) * (p_bounds - x0) / (x1 - x0)
}

get_refy_quantiles <- function(df, nobs = 2e4) {
  setorder(df, reporting_level, welfare)
  df_attr <- attributes(df)
  rls <- funique(df$reporting_level)
  p_bounds <- seq(0, 1, length.out = nobs + 1)

  qx <- lapply(rls, \(rl) {
    x <- df[reporting_level == rl & weight > 0]
    xpop <- fsum(x$weight)

    p_cum <- c(0, fcumsum(x$weight) / xpop)
    p_cum[length(p_cum)] <- 1
    cum_inc <- c(0, fcumsum(x$weight * x$welfare))

    interpolated_inc <- lorenz_interp(p_cum, cum_inc, p_bounds)

    bin_welfare <- diff(interpolated_inc) / (xpop / nobs)

    data.table(
      welfare = bin_welfare,
      weight = xpop / nobs,
      reporting_level = rl
    ) |>
      fselect(df_attr$names)
  }) |>
    rowbind()

  for (nm in names(df_attr)) {
    if (!nm %in% c("dim", "row.names", "names", "class", ".internal.selfref")) {
      setattr(qx, nm, df_attr[[nm]])
    }
  }

  qx
}
