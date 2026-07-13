/* / * / * / * / * / * / * / * / * / * / * / * / *
*
* Exercise 1: problematic countries and top-tail sensitivity
*
* This file uses pip_1kbins.dta, where mean is the PIP/survey mean
* and mean_1kbins is the mean reconstructed from the 1000-bin data.
*
* / * / * / * / * / * / * / * / * / * / * / * / */

/* Adjustable thresholds */
local rel_threshold = 0.01
local abs_threshold = 0.75

/* 1. Identify country-years with large differences */
use "$input\pip_1kbins.dta", clear

gen diff_mean = mean - mean_1kbins
gen abs_diff_mean = abs(diff_mean)
gen abs_diff_mean_pct = abs_diff_mean / mean

gen problematic = abs_diff_mean_pct >= `rel_threshold' | abs_diff_mean >= `abs_threshold'

gsort -abs_diff_mean_pct
order country_code country_name year welfare_type survey_acronym mean mean_1kbins ///
    diff_mean abs_diff_mean abs_diff_mean_pct problematic

export delimited country_code country_name year welfare_type survey_acronym ///
    mean mean_1kbins diff_mean abs_diff_mean abs_diff_mean_pct ///
    using "$path\work\problematic_countries.csv" if problematic, replace

preserve
    keep if problematic
    keep country_code country_name year welfare_type survey_acronym ///
        mean mean_1kbins diff_mean abs_diff_mean abs_diff_mean_pct pop
    rename pop pop_total
    tempfile problematic_cases
    save `problematic_cases', replace
restore

/* 2. Extract top-tail bins for those country-years */
use "$input\GlobalDist1000bins_1990_2026_20260324_2021_01_02_PROD.dta", clear

rename code country_code
merge m:1 country_code year using `problematic_cases', keep(match) nogen

rename pop pop_bin
keep if quantile >= 990

order country_code country_name year welfare_type survey_acronym quantile welf pop_bin ///
    mean mean_1kbins diff_mean abs_diff_mean abs_diff_mean_pct

export delimited using "$path\work\problematic_top_bins_990_1000.csv", replace

/* 3. Sensitivity: what top bin mean would exactly close the PIP vs 1k-bin gap? */
keep if quantile == 1000

gen topbin_share_of_1k_mean = (welf * pop_bin / pop_total) / mean_1kbins

/*
Because bin 1000 is about 0.1 percent of the population, changing only that bin
changes the total mean by:

    delta_total_mean = delta_top_bin_mean * pop_bin / pop_total

So the top-bin mean required to exactly match the PIP mean is:
*/
gen topbin_mean_if_match_survey = welf + (mean - mean_1kbins) * (pop_total / pop_bin)
gen topbin_pct_change_needed = topbin_mean_if_match_survey / welf - 1

order country_code country_name year welfare_type survey_acronym ///
    mean mean_1kbins diff_mean abs_diff_mean_pct welf pop_bin pop_total ///
    topbin_share_of_1k_mean topbin_mean_if_match_survey topbin_pct_change_needed

gsort -abs_diff_mean_pct

export delimited using "$path\work\problematic_topbin_sensitivity.csv", replace
