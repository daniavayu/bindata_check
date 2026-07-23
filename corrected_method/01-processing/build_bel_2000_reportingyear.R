library(haven)
library(data.table)

# Builds the processed data for BEL reporting year 2000.
# 2000 is itself a survey year (interpolated_means.dta has a single row for
# country_code == "BEL" & year == 2000, relative_distance == 1.0), so there
# is no interpolation/extrapolation between two surveys.
#
# Unlike MWI, BEL's raw data already comes pre-aggregated into 1,000
# population-weighted bins (one row per bin, not individual microdata).
# This script only does the PPP 2021 conversion + lineup mult_factor +
# weight rescale to the reporting population + bottom coding. The
# 1000 -> 20,000 bin expansion happens in 02-newmethod (each bin is
# duplicated 20x before applying the Lorenz-curve binning method).

ROOT  <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
INPUT <- file.path(ROOT, "corrected_method", "01-processing")
RAW   <- file.path(INPUT, "raw")

# 1. Lineup info for BEL, reporting year 2000 (single survey, no blending)
info <- as.data.table(read_dta(file.path(ROOT, "01-input", "interpolated_means.dta")))
info <- info[country_code == "BEL" & year == 2000 & welfare_time == 2000]
info[, year := 2000]
info <- info[, .(year, nac, nac_sy, relative_distance, predicted_mean_ppp, svy_mean, reporting_pop)]

cat("Lineup info (BEL, reporting year 2000):\n")
print(info)

# 2. Load raw 1,000-bin data
my_data <- as.data.table(read_dta(file.path(RAW, "BEL_2000.dta")))
my_data <- my_data[, .(bins, welfare, weight, cpi2021, icp2021, year)]

cat("\nRaw 1,000-bin data:\n")
cat("  Rows:", nrow(my_data), "\n")

# 3. Merge with lineup info by year
my_data <- merge(my_data, info, by = "year", all.x = TRUE)

# 4. Drop missing welfare
my_data <- my_data[!is.na(welfare)]

cat("\nAfter dropping missing welfare:\n")
cat("  Rows:", nrow(my_data), "\n")

# 5. Convert to 2021 PPP $/day
my_data[, welfare_ppp := welfare / (cpi2021 * icp2021 * 365)]

# 6. Lineup multiplicative factor (~1.0 here, since 2000 is its own survey year)
my_data[, mult_factor := predicted_mean_ppp / svy_mean]

# 7. Rescale population weights to the reporting population
my_data[, svy_pop := sum(weight)]
my_data[, weight := weight * (reporting_pop / svy_pop) * relative_distance]

# 8. Apply lineup factor to get final welfare
my_data[, welfare := welfare_ppp * mult_factor]

# 9. Bottom code at 0.28 (2021 PPP poverty-line floor)
my_data[welfare < 0.28, welfare := 0.28]

setorder(my_data, bins)

cat("\nFinal data (BEL, reporting year 2000, 1,000 bins):\n")
cat("  Rows:", nrow(my_data), "\n")
cat("  Mult factor:", unique(my_data$mult_factor), "\n")
cat("  Mean welfare:", weighted.mean(my_data$welfare, my_data$weight), "\n")

out_path <- file.path(INPUT, "bel_2000_reportingyear_processed.dta")
write_dta(my_data, out_path)
cat("\nSaved to:", out_path, "\n")
