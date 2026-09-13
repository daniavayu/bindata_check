# Scripts: what to run and what each file means

All scripts assume that R is started in the repository root (`bindata_check`).
General screening outputs go to `top_censoring/outputs/`; Pareto outputs go to
`top_censoring/pareto/outputs/`. Inputs are read from
`01-input/country/` and are never changed.

## Main analysis: Malawi 1997

| Run order | Script | Why it exists |
|---:|---|---|
| 1 | `01_lis_top_diagnostics.R` | Screens all surveys with the LIS log-IQR ceiling. |
| 2 | `../pareto/02_mwi_1997_pareto.R` | Visual upper-tail Pareto diagnostic, fitted from p90. |
| 3 | `../pareto/12_mwi_1997_eusilc_approaches.R` | **Main comparison table**: raw data, trimming, winsorizing, LIS, and robust-Pareto treatments. |

For the current working-paper task, script **12** is the most important file.
It uses survey weights both in the welfare statistics and in the Pareto fit.
It requires the R packages `haven` and `laeken`.

## Supporting Malawi diagnostics

| Script | Status | Purpose |
|---|---|---|
| `03_mwi_1997_index_sensitivity.R` | Supporting | Shows the diagnostic effect of excluding the two maximum observations. Not a proposed treatment. |
| `04_mwi_1997_dagum_sensitivity.R` | Supporting | Tests Dagum, a distribution for the full positive welfare distribution. |
| `05_mwi_1997_fisk_sensitivity.R` | Supporting | Tests Fisk/log-logistic, another full-distribution candidate. |
| `10_mwi_1997_treatment_scenarios.R` | Supporting | Earlier comparison of caps at the LIS ceiling. |
| `../pareto/11_mwi_1997_pareto_without_two.R` | Supporting | Confirms that the two maxima do not determine the broader Pareto-tail shape. |

These scripts explain *why* an adjustment may be needed, but are not necessary
to produce the main table.

## Screening and Belize analysis

| Script | Purpose |
|---|---|
| `06_screen_all_microdata.R` | Broad inventory of files and LIS diagnostics. |
| `07_prioritize_eligible_microdata.R` | Selects surveys with genuine microdata for this task. |
| `../pareto/08_blz_1993_1997_initial_analysis.R` | Initial Belize 1993 and 1997 tail analysis. |
| `../pareto/09_blz_pareto_threshold_sensitivity.R` | Checks whether p90, p95, or p99 is a better Pareto-tail starting point. |

## Naming rule for new scripts

Use a two-digit order plus country and purpose, for example:

`13_blz_1997_eusilc_approaches.R`

Create a new script only when it answers a new question. Otherwise, extend the
existing country comparison script.
