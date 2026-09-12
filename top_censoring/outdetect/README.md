# outdetect – Robust Top-Tail Outlier Screening (all GMD surveys)

## 1. Why this folder exists

The team already has a top-tail screening method based on the **LIS Tukey
fence** (see `../scripts/15_gmd_all_surveys_lis_indicator_comparison.py` and
`../outputs/24_surveys_four_indicators_lis_comparison.csv`), which flagged 24
survey-years as having a problematic concentration of welfare in a handful of
extreme top observations.

A separate request referencing the Stata
`outdetect` program (Vecchi/Mancini), used e.g. as:

```stata
outdetect qty_p if ..., norm(ln) zscore(median, s) alpha(3) out(top)
```

This folder implements the same *idea* — a robust, distribution-free outlier
identifier — in Python, applied to the same universe of surveys used by the
LIS method, so the two approaches can be compared apples-to-apples.

## 2. Method implemented

`outdetect` is a **robust z-score identifier**, confirmed against the source
paper: Belotti, F., Mancini, G. and Vecchi, G., *"Outlier Detection for
Welfare Analysis"*, World Bank Policy Research Working Paper
(https://openknowledge.worldbank.org/server/api/core/bitstreams/a8021d9e-bdcd-5371-b00c-589e10b5c0ec/content):

1. `norm(ln)` — work in `log(welfare)` space (reduces skewness before
   flagging).
2. `zscore(median, s)` — standardize each observation using a **robust**
   center and scale instead of the mean/SD (which are themselves distorted
   by outliers):
   - center = **unweighted median** of `log(welfare)`
   - scale = **unweighted S-estimator** (Rousseeuw & Croux, 1993) of
     `log(welfare)`: $S = 1.1926 \times \text{median}_i\{\text{median}_j |y_i - y_j|\}$.
     This is a *different* estimator from the MAD (Hampel, 1974) even though
     both are median-based — the paper's Table 2 lists `mad` and `s` as two
     separate scale options, and the requested syntax `zscore(median, s)`
     maps to the S row, not the MAD row.
3. `alpha(3)` — flag observations whose robust z-score exceeds `3`.
4. `out(top)` — only the **upper tail** is screened (this exercise is about
   top-censoring, not bottom outliers).

Threshold per survey:

$$
\text{threshold} = \exp\big(\text{median}(\ln w) + 3 \times 1.1926 \times \text{median}_i\{\text{median}_j |{\ln w_i - \ln w_j}|\}\big)
$$

Any record with `welfare > threshold` is flagged. Indicators
(`mean_ppp_day`, `gini`, `mld`, `atkinson_2`) are recomputed after capping
flagged values at the threshold ("outdetect-capped"), exactly mirroring how
the existing LIS-capped comparison table was built — these final indicators
**are weighted** (using each survey's sampling weights), for consistency
with the LIS-capped comparison.

### Center and scale are unweighted (confirmed from the paper)

Both the median (center) and the S-estimator (scale) are computed as
**unweighted** statistics, per the paper's explicit guidance:

> "we are assuming unweighted estimators. The weighted median does not, in
> fact, qualify as a robust statistic (Filzmoser, Gussenbauer and Templ
> 2016: 15)."

This is different from the LIS-ceiling method, whose 25th/75th percentiles
are computed **weighted**. Only the *final* indicators (`mean_ppp_day`,
`gini`, `mld`, `atkinson_2`, both "original" and "capped") are weighted for
both methods, for consistency with the sampling design.

### Computational note: subsampling for the S-estimator

The S-estimator is $O(n^2)$ — it requires the median of $n$ per-observation
medians of pairwise absolute differences, i.e. a full $n \times n$ pairwise
difference matrix. Survey sizes in this dataset range up to **5,020,899
records**, and 1,351 of the 1,825 surveys exceed 5,000 records, so an exact
computation is infeasible in Python (Stata's `outdetect` relies on a fast
Mata implementation of the Croux–Rousseeuw (1992) algorithm, which is not
replicated here).

As a practical compromise, when a survey has more than
`S_ESTIMATOR_SUBSAMPLE_SIZE = 2000` records, the S-estimator (only) is
computed on a fixed-seed (`RANDOM_SEED = 42`) random subsample of 2,000
records drawn from that survey. The median center, and the resulting
threshold/capping/indicator computations, are still applied to the **full**
survey. The output column `S-estimator subsampled (Y/N)` in
`outputs/outdetect_all_surveys_summary.csv` flags which surveys used this
compromise (1,351 of 1,825).

### Grouping

Exactly like the LIS method: **one screening per full survey (Country +
Year + Survey)**, not per subgroup — confirmed with the requester as the
right level of granularity to match the existing comparison.

## 3. Universe covered

Both scripts read the full `01-input/GMD_all_2017.dta` file and apply the
same microdata filter used everywhere else in this repo:

- `use_microdata == 1`
- `use_bin == 0 & use_groupdata == 0` (i.e., not aggregated/grouped data)
- `welfare > 0` and finite
- `weight > 0` and finite

This produces the same **1,825 survey-year rows** as
`../outputs/all_surveys_four_indicators_lis_comparison.csv`, so the two
methods are directly comparable.

## 4. Scripts

| Script | Purpose | Output |
|---|---|---|
| [`01_outdetect_all_surveys.py`](01_outdetect_all_surveys.py) | Computes the outdetect (median+S-estimator, log, alpha=3, top-only) threshold and capped indicators for all 1,825 surveys. | `outputs/outdetect_all_surveys_summary.csv` |
| [`02_compare_lis_vs_outdetect.py`](02_compare_lis_vs_outdetect.py) | Merges the outdetect results with the existing LIS results and flags surveys caught by each method. | `outputs/lis_vs_outdetect_comparison.csv`, `outputs/always_problematic_surveys.csv` |
| [`03_executive_summary.py`](03_executive_summary.py) | Recomputes weighted indicators for the current intersection and creates the Malawi 1997 diagnostic table. | `outputs/executive_summary_always_problematic.csv`, `outputs/executive_summary_mwi_1997_flagship_top_values.csv` |

Run order: `01_outdetect_all_surveys.py`, then `02_compare_lis_vs_outdetect.py`,
then `03_executive_summary.py`
(from the `bindata_check` directory, using the repo's Python environment).

## 5. How "problematic" is defined for the comparison

- **LIS-flagged** = survey appears in the existing, team-confirmed
  `24_surveys_four_indicators_lis_comparison.csv`. Verified rule: exactly the
  surveys with `Welfare share above ceiling (%) > 5%` in the full universe
  (min value among the 24 is PRY 2010 at 5.09%; the next survey down, UGA
  2019, is at 4.86% and is excluded) — a clean 5% cutoff, confirmed by
  re-deriving it directly from `all_surveys_four_indicators_lis_comparison.csv`.
- **outdetect-flagged** = at least **5% of total welfare** is held by
  observations flagged by outdetect. This is the same welfare-share criterion
  used for the LIS screening and flags **398 of 1,825 surveys**.
  The 5% threshold is an empirical comparability criterion adopted for this
  exercise; it is not an official survey-level cutoff prescribed by the
  outdetect documentation.
- **Flagged by both ("always problematic")** = the survey shows up as an
  extreme top-tail case *regardless of which of the two outlier-screening
  methods is used* — the strongest possible evidence that these surveys need
  attention.

## 6. Results

Of the 24 LIS-flagged surveys and the 398 outdetect-flagged surveys, drawn
from the same 1,825-survey universe, **24 surveys are flagged by both
methods**. These are the surveys where an extreme top tail is identified by
both the Tukey-IQR fence and the S-estimator (Rousseeuw & Croux, 1993):

| Country | Year | Survey | Welfare type | LIS welfare share above ceiling (%) | outdetect welfare share flagged (%) |
|---|---|---|---|---|---|
| MWI | 1997 | IHS-I | consumption | 41.60 | 47.27 |
| BLZ | 1993 | LFS | income | 20.09 | 21.49 |
| BLZ | 1994 | LFS | income | 19.96 | 19.96 |
| BLZ | 1997 | LFS | income | 18.24 | 19.43 |
| SYC | 2013 | HBS | income | 12.06 | 17.61 |
| MWI | 2016 | IHS-IV | consumption | 10.97 | 15.63 |
| ETH | 1995 | HICES | consumption | 10.89 | 15.74 |
| COM | 2004 | EIM | consumption | 7.61 | 22.30 |
| RWA | 2005 | EICV-II | consumption | 6.14 | 16.27 |
| RWA | 2013 | EICV-IV | consumption | 5.37 | 14.85 |

No surveys are flagged only by LIS under this criterion. The 374 surveys
flagged only by outdetect are listed in
`outputs/lis_vs_outdetect_comparison.csv` (columns `Flagged by LIS`,
`Flagged by outdetect`) for the complete breakdown.

**Takeaway:** all 24 LIS-flagged surveys also exceed the 5% welfare-share
threshold under outdetect. The outdetect criterion is substantially more
sensitive and identifies 374 additional surveys, so the two methods agree on
the LIS cases but outdetect produces a much broader screening list.

See `outputs/always_problematic_surveys.csv` for the full list with values,
and `outputs/lis_vs_outdetect_comparison.csv` for the full 1,825-row
side-by-side comparison (includes `welfare_type`, consumption vs. income,
merged from `interpolated_means.csv`).
