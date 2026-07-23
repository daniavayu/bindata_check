library(haven)
library(data.table)

# Builds the pooled microdata for MWI reporting year 2004.
# Reporting year 2004 falls BETWEEN two surveys (1997 and 2004), so both are
# pooled together, each reweighted by how close it is to 2004
# (relative_distance), per interpolated_means.dta.

ROOT <- "C:/Users/wb661551/OneDrive - WBG/Desktop/Internship/Bottom Censoring/bindata_check"
INPUT <- file.path(ROOT, "corrected_method", "01-processing")
RAW   <- file.path(INPUT, "raw")

# 1. Lineup info for MWI, reporting year 2004: blends survey years 1997.83 & 2004.23
info <- as.data.table(read_dta(file.path(ROOT, "01-input", "interpolated_means.dta")))
info <- info[country_code == "MWI" & year == 2004 & welfare_time %in% c(1997.83, 2004.23)]
info[, year := fifelse(welfare_time == 1997.83, 1997, 2004)]
info <- info[, .(year, nac, nac_sy, relative_distance, predicted_mean_ppp, svy_mean, reporting_pop)]

cat("Lineup info (reporting year 2004):\n")
print(info)

# 2. Load raw microdata and append (2004 first, then 1997 - matches .do order)
d2004 <- as.data.table(read_dta(file.path(RAW, "MWI_2004.dta")))
setnames(d2004, "weight_p", "weight")

d1997 <- as.data.table(read_dta(file.path(RAW, "MWI_1997.dta")))

my_data <- rbindlist(list(d2004, d1997), fill = TRUE)
my_data <- my_data[, .(welfare, weight, cpi2021, icp2021, year)]

cat("\nAfter append:\n")
cat("  Rows:", nrow(my_data), "\n")

# 3. Merge with lineup info by year
my_data <- merge(my_data, info, by = "year", all.x = TRUE)

# 4. Drop missing welfare
my_data <- my_data[!is.na(welfare)]

cat("\nAfter dropping missing welfare:\n")
cat("  Rows:", nrow(my_data), "\n")

# 5. Convert to 2021 PPP $/day
my_data[, welfare_ppp := welfare / (cpi2021 * icp2021 * 365)]

# 6. Lineup multiplicative factor
my_data[, mult_factor := predicted_mean_ppp / svy_mean]

# 7. Rescale population weights (by year) to the 2004 reporting population,
#    weighted by relative_distance (1997 survey contributes little, ~3.6%)
my_data[, svy_pop := sum(weight), by = year]
my_data[, weight := weight * (reporting_pop / svy_pop) * relative_distance]

# 8. Apply lineup factor to get final welfare
my_data[, welfare := welfare_ppp * mult_factor]

# 9. Bottom code at 0.28 (2021 PPP poverty-line floor)
my_data[welfare < 0.28, welfare := 0.28]

cat("\nFinal pooled data (reporting year 2004):\n")
cat("  Rows:", nrow(my_data), "\n")
cat("  Mean welfare:", weighted.mean(my_data$welfare, my_data$weight), "\n")

out_path <- file.path(INPUT, "mwi_2004_reportingyear_processed.dta")
write_dta(my_data, out_path)
cat("\nSaved to:", out_path, "\n")
