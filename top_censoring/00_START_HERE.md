# Top-censoring: start here

This folder contains the separate analysis of whether very high welfare values
should be treated before reporting inequality measures. 

## The one-sentence question

Do a few observations in the upper tail have enough influence on mean welfare
and inequality that a transparent top-coding treatment is warranted?

## Recommended route: Malawi 1997

You do **not** need to run every script. Start with these three:

1. `scripts/01_lis_top_diagnostics.R`
   - Applies the LIS log-IQR rule to screen all available surveys.
2. `scripts/02_mwi_1997_pareto.R`
   - Makes the Pareto tail graph for Malawi 1997.
3. `scripts/12_mwi_1997_eusilc_approaches.R`
   - Produces the main comparison of no adjustment, trimming, winsorizing,
     LIS top-coding, and Pareto-based treatments.

Run from the repository root:

```r
source("top_censoring/scripts/12_mwi_1997_eusilc_approaches.R")
```

The main table is:

`outputs/mwi_1997_eusilc_approaches_table.csv`

## Current conclusion for Malawi 1997

- A small upper-tail group strongly changes the mean and inequality measures.
- Poverty at the $2.41/day 2021-PPP line is essentially unchanged.
- LIS and Pareto-based treatments produce similar conclusions: some treatment
  of the upper tail is warranted for sensitivity and reporting purposes.
- The evidence supports **top-coding, not deleting observations**.

## Folder map

| Folder or file | Use it for |
|---|---|
| `scripts/` | Reproducible R analysis. See its README before running files. |
| `notebooks/` | Interactive exploration of one survey at a time. |
| `outputs/` | Generated tables and figures. See its README for the important files. |
| `RESULTS_MWI_1997.md` | Earlier Malawi findings and interpretation. |
| `WORKFLOW.md` | Original workflow notes; use this file and `scripts/README.md` as the current guide. |

## What not to do

- Do not fit Pareto to the entire welfare distribution: it models only the
  upper tail.
- Do not delete extreme records in production estimates just because they are
  high. Deletion is useful only as a diagnostic sensitivity check.
- Do not use binned data for this exercise. Top-outlier identification requires
  the original survey microdata and its person weights.
