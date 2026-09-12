# Outputs: what to look at

Outputs are generated files. They can be recreated by running the scripts in
`../scripts/` and `../pareto/`; do not edit them manually.

Pareto outputs are stored in the `pareto/` subfolder.

## Main Malawi 1997 outputs

| File | Read this when you need |
|---|---|
| `pareto/mwi_1997_eusilc_approaches_table.csv` | One row per treatment and one column per indicator. This is the main comparison table. |
| `pareto/mwi_1997_eusilc_approaches_comparison.csv` | The same results in long form, including percentage change from raw data. |
| `pareto/mwi_1997_eusilc_method_notes.csv` | Plain-language description of every approach. |
| `pareto/mwi_1997_eusilc_settings.csv` | Exact thresholds, Pareto parameter, and number of flagged records. |
| `pareto/mwi_1997_eusilc_imputation_uncertainty.csv` | Simulation variation from Pareto multiple imputation. |
| `pareto/mwi_1997_pareto_p90_ccdf.png` | The first Pareto diagnostic graph. |

## Other groups of outputs

- Files beginning `blz_`: Belize 1993 and 1997 analyses.
- Files beginning `all_microdata`, `eligible_microdata`, `survey_inventory`, or
  `lis_top`: survey-screening results.
- Files beginning `mwi_1997_dagum`, `mwi_1997_fisk`, or
  `mwi_1997_pareto_without`: supporting distribution-fit exercises.
- Files beginning `mwi_1997_two_extremes` or `mwi_1997_treatment_scenarios`:
  earlier diagnostic sensitivity exercises.
