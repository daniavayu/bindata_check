# Pareto analyses

This folder contains upper-tail analyses based on Pareto models.

## Scripts

- `02_mwi_1997_pareto.R`: Malawi 1997 p90 Pareto diagnostic.
- `08_blz_1993_1997_initial_analysis.R`: Belize 1993 and 1997 Pareto diagnostics and index sensitivity.
- `09_blz_pareto_threshold_sensitivity.R`: Belize threshold sensitivity.
- `11_mwi_1997_pareto_without_two.R`: Malawi Pareto sensitivity excluding the two largest observations.
- `12_mwi_1997_eusilc_approaches.R`: Malawi comparison of trimming, winsorizing, and robust Pareto treatments.

Run from the `bindata_check` directory. Outputs are written to `../outputs/pareto/`.

Pareto is one of three conceptually separate approaches in this project:

1. Pareto: whether the observed upper tail follows a plausible fitted tail.
2. GB2: whether the full distribution or body-implied tail can reproduce the extremes.
3. Regression: whether observed extremes can be explained by household characteristics.
