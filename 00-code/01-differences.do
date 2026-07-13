/* / * / * / * / * / * / * / * / * / * / * / * / * 
*
* 2. Identify large differences
*
* / * / * / * / * / * / * / * / * / * / * / * / */

use "$input\pip_1kbins.dta", clear

gen diff_mean = abs(mean_1kbins - mean)
gsort -diff_mean

br country_code year diff_mean mean mean_1kbins

gen diff_mean_pct = diff_mean/mean
gsort -diff_mean_pct
br country_code year diff_mean_pct mean mean_1kbins
