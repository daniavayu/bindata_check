# Top-censoring diagnostics

This directory evaluates whether a consistent top-coding rule is warranted for
the available *observed survey microdata*. This directory contains the
top-censoring analysis, separate from the bottom-coding and bin-reconstruction
analyses.

**Start here:** [`00_START_HERE.md`](00_START_HERE.md) for the recommended
route through the Malawi 1997 case study. This file documents the overall
scope and folder map.

## Scope

`01-input/country` currently contains two input types:

- **Microdata:** `AGO_2008.dta`, `MWI_1997.dta`, `MWI_2004.dta`, and
  `MWI_2019.dta`. These are eligible for the LIS outlier and distribution-fit
  analyses.
- **Grouped distributions:** `BEL_2000.dta` and `CHE_2000.dta`, identified by
  their `bins` column. They may be used only for descriptive sensitivity checks;
  grouped values cannot establish that individual observations are top outliers.

The analysis also draws on `01-input/GMD_all_2017.dta`, the full ~1,825-survey
GMD universe, for cross-survey screening (LIS and outdetect methods).

The analysis uses each survey's welfare aggregate and its person weight. Welfare
is also expressed as 2021 PPP USD per person per day using
`welfare / (cpi2021 * icp2021 * 365)` when the conversion fields are available.
The LIS classification is invariant to this within-survey rescaling.

## Analyses in this folder

| Folder | Method | Details |
|---|---|---|
| `scripts/` | LIS log-IQR screening (single surveys and full GMD universe) | [`scripts/README.md`](scripts/README.md) |
| `pareto/` | Pareto/Dagum/Fisk upper-tail and full-distribution fits | folder README |
| `gb2/` | GB2 distribution fit, full sample vs. top observations removed | [`gb2/README.md`](gb2/README.md) |
| `regression/` | Predicted-vs-observed welfare using household characteristics, cross-fitted residual flags | [`regression/README.md`](regression/README.md) |
| `outdetect/` | Robust top-tail screening (median + S-estimator z-score, Belotti/Mancini/Vecchi `outdetect`), applied to all 1,825 GMD surveys and compared against the LIS method. S-estimator subsampled for large surveys (Python, O(n^2)). | [`outdetect/README.md`](outdetect/README.md) |
| `outdetect_r/` | Same outdetect method, but with an **exact** (non-subsampled) S-estimator via R's `robustbase::Sn()` (O(n log n), same algorithm Stata's `outdetect` uses internally) — no random seed, no approximation, for all 1,825 surveys. | [`outdetect_r/README.md`](outdetect_r/README.md) |

## Sequence

1. Run `scripts/01_lis_top_diagnostics.R` to inventory the files and calculate
   the weighted LIS ceiling on log welfare.
2. Inspect the descriptive results. A small count alone is not evidence that
   top-coding is necessary; examine weighted population share, welfare share,
   maximum-to-ceiling ratio, and sensitivity of inequality measures.
3. Fit candidate distributions to the original, uncensored microdata:
   Pareto for an explicitly selected upper-tail threshold, GB2 for a
   flexible four-parameter fit, and Dagum/Fisk for the full positive
   distribution.
4. Cross-check with `outdetect/` (an independent, robust-statistics screening
   method) and `regression/` (an independent, covariate-based screening
   method) to see whether the same observations are flagged by unrelated
   approaches.
5. Treat LIS top-coding as a sensitivity analysis. Refit after top-coding only
   to describe robustness, never as the basis for deciding that the removed
   observations were outliers.

## LIS ceiling

For strictly positive welfare `y`, calculate weighted quartiles of `log(y)`:

`ceiling = exp(Q3(log(y)) + k * (Q3(log(y)) - Q1(log(y))))`.

LIS uses `k = 3`. The script also reports `k = 2` and `k = 4` as sensitivity
benchmarks. Observations above the ceiling are candidates for review, not
automatically invalid data.

## Outputs

`scripts/01_lis_top_diagnostics.R` writes non-versioned CSVs to
`top_censoring/outputs/`:

- `survey_inventory.csv`
- `lis_top_diagnostics.csv`

Each of the other analyses writes to its own `outputs/` subfolder
(`pareto/outputs/`, `gb2/outputs/`, `regression/outputs/`, `outdetect/outputs/`)
— see the corresponding folder README for the specific files and how to
interpret them.

No source microdata are modified by any script in this folder.
