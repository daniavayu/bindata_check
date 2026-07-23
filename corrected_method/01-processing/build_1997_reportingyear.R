library(haven)
library(data.table)

# Builds the microdata for MWI reporting year 1997.
# 1997 is itself a survey year (interpolated_means.dta has a single row for
# year == 1997, relative_distance == 1.0), so there is NO interpolation or
# extrapolation between two surveys - unlike reporting year 2004, which
# blends the 1997 and 2004 surveys.
#
# Still goes through: PPP 2021 conversion + lineup mult_factor + weight
# rescale to the 1997 reporting population + bottom coding at 0.28.

ROOT  <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
INPUT <- file.path(ROOT, "corrected_method", "01-processing")
RAW   <- file.path(INPUT, "raw")

# 1. Lineup info for MWI, reporting year 1997 (single survey, no blending)
info <- as.data.table(read_dta(file.path(ROOT, "01-input", "interpolated_means.dta")))
info <- info[country_code == "MWI" & year == 1997 & welfare_time == 1997.83]
info[, year := 1997]
info <- info[, .(year, nac, nac_sy, relative_distance, predicted_mean_ppp, svy_mean, reporting_pop)]

cat("Lineup info (reporting year 1997):\n")
print(info)

# 2. Load raw microdata for the 1997 survey only
my_data <- as.data.table(read_dta(file.path(RAW, "MWI_1997.dta")))
my_data <- my_data[, .(welfare, weight, cpi2021, icp2021, year)]

cat("\nRaw microdata:\n")
cat("  Rows:", nrow(my_data), "\n")

# 3. Merge with lineup info by year
my_data <- merge(my_data, info, by = "year", all.x = TRUE)

# 4. Drop missing welfare
my_data <- my_data[!is.na(welfare)]

cat("\nAfter dropping missing welfare:\n")
cat("  Rows:", nrow(my_data), "\n")

# 5. Convert to 2021 PPP $/day
my_data[, welfare_ppp := welfare / (cpi2021 * icp2021 * 365)]

# 6. Lineup multiplicative factor (predicted mean for reporting year 1997 /
#    the 1997 survey's own mean - no blending with another survey)
my_data[, mult_factor := predicted_mean_ppp / svy_mean]

# 7. Rescale population weights to the 1997 reporting population
#    (relative_distance == 1, so this is effectively weight * reporting_pop/svy_pop)
my_data[, svy_pop := sum(weight)]
my_data[, weight := weight * (reporting_pop / svy_pop) * relative_distance]

# 8. Apply lineup factor to get final welfare
my_data[, welfare := welfare_ppp * mult_factor]

# 9. Bottom code at 0.28 (2021 PPP poverty-line floor)
my_data[welfare < 0.28, welfare := 0.28]

cat("\nFinal data (reporting year 1997):\n")
cat("  Rows:", nrow(my_data), "\n")
cat("  Mult factor:", unique(my_data$mult_factor), "\n")
cat("  Mean welfare:", weighted.mean(my_data$welfare, my_data$weight), "\n")

out_path <- file.path(INPUT, "mwi_1997_reportingyear_processed.dta")
write_dta(my_data, out_path)
cat("\nSaved to:", out_path, "\n")
