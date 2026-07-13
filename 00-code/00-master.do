/* / * / * / * / * / * / * / * / * / * / * / * / * 
*
* 1. Set paths and save merged data
*
* / * / * / * / * / * / * / * / * / * / * / * / */

/* Set paths */
gl path "C:\Users\wb661551\OneDrive - WBG\Desktop\Internship\Bottom Censoring"

gl input "$path\01-input"

*1. pip data
pip, clear fillgaps 

    keep country_code year welfare_type reporting_level mean
    tempfile all
    save `all', replace

pip, clear 
    drop mean
    merge 1:1 country_code year welfare_type reporting_level using `all', keep(3) nogen

    drop if reporting_level!="national" & country_code=="CHN"

save "$input\1-pip.dta", replace


*2. Create means from 1k bin data
use "$input\GlobalDist1000bins_1990_2026_20260324_2021_01_02_PROD.dta", clear

//headcount 
foreach pov in 3 /*4.2 8.3*/ {
    local x=`pov'*10
    gen hc_`x' = welf<`pov'
}

sort code year quantile

collapse (mean) welf hc* (rawsum) pop [w=pop], by(code year)
ren hc_ headcount_1kbins
ren welf mean_1kbins

sort code year 
ren code country_code 

tempfile data_1kbins
save `data_1kbins', replace


*3. merge data
use "$input\1-pip.dta", clear
merge 1:1 country_code year using `data_1kbins', keep(3) nogen

save "$input\pip_1kbins.dta", replace

do "$path\00-code\02-problematic-top-tail.do"
