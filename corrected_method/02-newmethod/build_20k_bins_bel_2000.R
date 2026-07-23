library(haven)
library(data.table)
library(fastverse)
library(openxlsx)

fastverse_conflicts()
conflicted::conflicts_prefer(data.table::between())
conflicted::conflicts_prefer(collapse::funique())
conflicted::conflicts_prefer(data.table::first())

# New method: builds 20,000 bins for BEL reporting year 2000.
#
# Unlike MWI (individual microdata), BEL's processed data is already
# 1,000 population-weighted bins. There is no extra information about the
# distribution *within* each bin, so to go from 1,000 -> 20,000 bins we
# duplicate each bin 20x (same welfare value, weight split into 20 equal
# parts). This preserves the exact mean and the exact shape of the
# distribution (same Lorenz curve), just at a finer (20k) resolution.
#
# We then run the same generalized-Lorenz-curve binning function used for
# MWI on this expanded pseudo-microdata, for full consistency of method
# and as a validation step (it should return ~identical bin values).

ROOT   <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
INPUT  <- file.path(ROOT, "corrected_method", "01-processing")
OUTPUT <- file.path(ROOT, "corrected_method", "03-outputs")
bel_path <- file.path(INPUT, "bel_2000_reportingyear_processed.dta")

# Read the 1,000-bin processed data
bel_1000 <- as.data.table(read_dta(bel_path))

cat("Input data (BEL, reporting year 2000, 1,000 bins):\n")
cat("  Rows:", nrow(bel_1000), "\n")
cat("  Weighted mean:", weighted.mean(bel_1000$welfare, bel_1000$weight), "\n\n")

# 1. Expand: duplicate each bin 20x, splitting its weight into 20 equal parts
n_copies <- 20
my_data <- bel_1000[rep(1:.N, each = n_copies)]
my_data[, weight := weight / n_copies]
my_data[, reporting_level := "national"]

cat("After expanding to pseudo-microdata (1,000 x 20):\n")
cat("  Rows:", nrow(my_data), "\n")
cat("  Weighted mean:", weighted.mean(my_data$welfare, my_data$weight), "\n\n")

get_refy_quantiles_mean <- function(df, nobs = 2e4) {

  setorder(df, welfare)
  df_attr <- attributes(df)
  rls <- funique(df$reporting_level)

  # 1. Define the 20,001 exact boundaries of the bins (from 0 to 1)
  p_bounds <- seq(0, 1, length.out = nobs + 1)

  qx <- lapply(rls, \(rl) {
    # Filter out 0-weight rows to avoid breaking the cumulative interpolation
    x <- df[reporting_level == rl & weight > 0]
    xpop <- fsum(x$weight)

    # 2. Calculate cumulative population share (p_cum)
    p_cum <- c(0, fcumsum(x$weight) / xpop)
    p_cum[length(p_cum)] <- 1  # guard against floating-point rounding (e.g. 0.999999...)

    # 3. Calculate cumulative total welfare (Generalized Lorenz curve)
    cum_inc <- c(0, fcumsum(x$weight * x$welfare))

    # 4. Interpolate the exact cumulative welfare at the 20,000 boundaries
    # Linear interpolation perfectly splits households that overlap two bins
    interpolated_inc <- approx(x = p_cum, y = cum_inc, xout = p_bounds, method = "linear")$y

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

cat("Creating 20k bins...\n")
result <- get_refy_quantiles_mean(my_data, nobs = 2e4)

cat("\nResult:\n")
cat("  Rows:", nrow(result), "\n")
cat("  Weighted mean:", weighted.mean(result$welfare, result$weight), "\n")
cat("  Total weight:", sum(result$weight), "\n")

output_path <- file.path(OUTPUT, "BEL_2000_reportingyear_20k_mean.xlsx")
write.xlsx(result, output_path)
cat("\nSaved to:", output_path, "\n")
