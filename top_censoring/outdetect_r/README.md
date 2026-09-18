# outdetect_r – Exact (non-subsampled) outdetect screening

## Why this folder exists

[`../outdetect/`](../outdetect/) implements the `outdetect` robust top-tail
screening method in Python, but the S-estimator (Rousseeuw & Croux, 1993)
is O(n^2), so it has to be computed on a fixed-seed random **subsample**
of at most 2,000 records for any survey larger than that (1,351 of the
1,825 GMD surveys). That is a reasonable practical compromise, but it means
the reported number for a given survey is not fully independent of every
other survey processed before it in the same run (they share one random
number generator). This shared random-number stream previously affected the
Malawi 1997 results because the subsample could depend on the processing order
of the surveys.

This folder removes the compromise entirely. It replicates the exact same
`outdetect welfare, norm(ln) zscore(median, s) alpha(3) out(top)` method,
but computes the S-estimator **exactly**, for every survey, using R's
`robustbase::Sn()` — an implementation of the Croux & Rousseeuw (1992)
O(n log n) algorithm, the same fast algorithm Stata's `outdetect` uses
internally in Mata. No subsampling, no random seed, anywhere in this folder.

## Why this is feasible (benchmark)

| n | naive O(n^2) (Python, this repo) | `robustbase::Sn()` (O(n log n)) |
|---|---|---|
| 2,000 | ~0.2 s | ~0.001 s |
| 20,000 | ~35 s (chunked) | ~0.009 s |
| 200,000 | infeasible in practice | ~0.11 s |
| 1,000,000 | infeasible in practice | ~0.49 s |

Even the largest GMD survey (up to ~5 million records) is computed exactly
in well under a second. The full ~1,825-survey screening finishes in a
couple of minutes.

## Correctness check

`Sn(x, constant = 1.1926)` was checked against this repo's naive formula
`1.1926 * median_i(median_j|x_i - x_j|)` on random samples:

| n | naive formula | `Sn()` |
|---|---|---|
| 50 | 1.166362 | 1.188843 |
| 200 | 0.921066 | 0.924732 |
| 1,000 | 0.987137 | 0.987427 |

The (shrinking, as n grows) difference is expected: `Sn()` implements the
precise order-statistic definition of Sn from Rousseeuw & Croux (1993),
which is more faithful to the original estimator than the naive nested-median
approximation — not a discrepancy to be concerned about.

## Requirements

```r
install.packages("robustbase")
```

(installed in this environment already; see the correctness/speed checks
above for how it was validated before adopting it here).

## Scripts

Run from the `bindata_check` directory, in this order:

| Script | Purpose | Output |
|---|---|---|
| [`01_outdetect_all_surveys.R`](01_outdetect_all_surveys.R) | Exact outdetect threshold and capped indicators for all ~1,825 GMD surveys. | `outputs/outdetect_all_surveys_summary.csv` |
| [`02_compare_lis_vs_outdetect.R`](02_compare_lis_vs_outdetect.R) | Merges with the existing LIS results; flags surveys caught by each method. | `outputs/lis_vs_outdetect_comparison.csv`, `outputs/always_problematic_surveys.csv` |
| [`03_executive_summary.R`](03_executive_summary.R) | Recomputes weighted indicators for the "always problematic" surveys and builds the Malawi 1997 deep-dive tables (top-30 values, 5-scenario treatment comparison). | `outputs/executive_summary_always_problematic.csv`, `outputs/executive_summary_mwi_1997_flagship_top_values.csv`, `outputs/mwi_1997_treatment_comparison_5_scenarios.csv` |

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" top_censoring\outdetect_r\01_outdetect_all_surveys.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" top_censoring\outdetect_r\02_compare_lis_vs_outdetect.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" top_censoring\outdetect_r\03_executive_summary.R
```

## How this differs from `../outdetect/`

- Method, thresholds, flagging rule (welfare share >= 5%), and output schema
  are the same as the Python version — see [`../outdetect/README.md`](../outdetect/README.md)
  for the full methodology writeup (unweighted median/S-estimator in log
  space, weighted final indicators, alpha = 3, top-tail only).
- The only substantive difference is **how** the S-estimator is computed:
  exact here (`robustbase::Sn()`), subsampled there (naive O(n^2) on a
  2,000-record random draw for large surveys). Because of this, the exact
  numbers in `outputs/` here will differ slightly from `../outdetect/outputs/`
  for large surveys, and will match closely (not identically, since `Sn()`
  uses the precise order-statistic definition, see above) for small ones.
- This version has no `RANDOM_SEED` or `S_ESTIMATOR_SUBSAMPLE_SIZE` constant,
  and no "S-estimator subsampled (Y/N)" column, because none of that is
  needed anymore.
