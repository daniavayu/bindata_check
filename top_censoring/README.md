# Top-censoring diagnostics

This directory evaluates whether a consistent top-coding rule is warranted for
the available *observed survey microdata*. It is deliberately separate from the
completed bottom-coding and bin-reconstruction work.

## Scope

`01-input/country` currently contains two input types:

- **Microdata:** `AGO_2008.dta`, `MWI_1997.dta`, `MWI_2004.dta`, and
  `MWI_2019.dta`. These are eligible for the LIS outlier and distribution-fit
  analyses.
- **Grouped distributions:** `BEL_2000.dta` and `CHE_2000.dta`, identified by
  their `bins` column. They may be used only for descriptive sensitivity checks;
  grouped values cannot establish that individual observations are top outliers.

The analysis uses each survey's welfare aggregate and its person weight. Welfare
is also expressed as 2021 PPP USD per person per day using
`welfare / (cpi2021 * icp2021 * 365)` when the conversion fields are available.
The LIS classification is invariant to this within-survey rescaling.

## Sequence

1. Run `scripts/01_lis_top_diagnostics.R` to inventory the files and calculate
   the weighted LIS ceiling on log welfare.
2. Inspect the descriptive results. A small count alone is not evidence that
   top-coding is necessary; examine weighted population share, welfare share,
   maximum-to-ceiling ratio, and sensitivity of inequality measures.
3. Fit candidate distributions to the original, uncensored microdata:
   Pareto for an explicitly selected upper-tail threshold and Dagum/Fisk for
   the full positive distribution.
4. Treat LIS top-coding as a sensitivity analysis. Refit after top-coding only
   to describe robustness, never as the basis for deciding that the removed
   observations were outliers.

## LIS ceiling

For strictly positive welfare `y`, calculate weighted quartiles of `log(y)`:

`ceiling = exp(Q3(log(y)) + k * (Q3(log(y)) - Q1(log(y))))`.

LIS uses `k = 3`. The script also reports `k = 2` and `k = 4` as sensitivity
benchmarks. Observations above the ceiling are candidates for review, not
automatically invalid data.

## Outputs

The script writes non-versioned CSVs to `top_censoring/outputs/`:

- `survey_inventory.csv`
- `lis_top_diagnostics.csv`

No source microdata are modified.
