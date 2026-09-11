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

Run order: `01_outdetect_all_surveys.py` then `02_compare_lis_vs_outdetect.py`
(from the `bindata_check` directory, using the repo's Python environment).

## 5. How "problematic" is defined for the comparison

- **LIS-flagged** = survey appears in the existing, team-confirmed
  `24_surveys_four_indicators_lis_comparison.csv`. Verified rule: exactly the
  surveys with `Welfare share above ceiling (%) > 5%` in the full universe
  (min value among the 24 is PRY 2010 at 5.09%; the next survey down, UGA
  2019, is at 4.86% and is excluded) — a clean 5% cutoff, confirmed by
  re-deriving it directly from `all_surveys_four_indicators_lis_comparison.csv`.
- **outdetect-flagged** = survey is in the **top 24** surveys (by "Welfare
  share flagged (%)") under the outdetect method (26 survey-rows pass this
  cut in practice, due to ties at the 24th-place value).
  **A flat 5% cutoff does *not* carry over to outdetect**: the two methods
  flag very different baseline shares of welfare per survey
  (LIS mean/median across all 1,825 surveys = 0.85% / 0.38%; outdetect
  mean/median = 3.57% / 2.80%, since the S-estimator identifier is inherently
  more sensitive than the Tukey/IQR fence). Applying LIS's raw 5% threshold to
  outdetect's "Welfare share flagged (%)" column flags **398 of 1,825
  surveys** (~22%), which is not a usable "problematic" list.
  Instead, the LIS 5% cutoff was translated to its **equivalent percentile
  rank** (24/1825 flagged → the 98.685th percentile of the LIS distribution),
  and that same percentile rank was applied to the outdetect distribution.
  That yields an outdetect threshold of **`Welfare share flagged (%) > 13.77%`**,
  close to taking the top 24 (26, with ties) by that column, which is what
  `02_compare_lis_vs_outdetect.py` does directly (`N_TOP = 24`). So the "24"
  is not a coincidence and not an arbitrary tie-in to the LIS list size
  either — it is what falls out of applying an equally-strict (same-percentile)
  cutoff to both methods.
- **Flagged by both ("always problematic")** = the survey shows up as an
  extreme top-tail case *regardless of which of the two outlier-screening
  methods is used* — the strongest possible evidence that these surveys need
  attention.

## 6. Results

Of the 24 LIS-flagged surveys and the 26 outdetect-flagged surveys (ties at
the 24th-place cutoff), drawn from the same 1,825-survey universe,
**10 surveys are flagged by both methods** — these are the surveys where an
extreme top tail shows up regardless of whether a Tukey-IQR fence or an
S-estimator (Rousseeuw & Croux, 1993) identifier is used:

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

14 surveys are flagged only by LIS (Tukey/IQR) and not by outdetect:
MDG 1993, MDV 2002, MMR 2015, NOR 2005, NOR 2006, PRY 2002, PRY 2007,
PRY 2010, PRY 2014, MOZ 2014, BFA 1998, MOZ 2002, CAF 2008, RWA 2000.
16 different surveys are flagged only by outdetect (S-estimator) and not by
LIS, notably a cluster of Russia HBS/VNDN survey-years (2007–2019) and Chile
CASEN survey-years (1996, 2000, 2009) — see
`outputs/lis_vs_outdetect_comparison.csv` (columns `Flagged by LIS`,
`Flagged by outdetect`) for the complete breakdown.

**Takeaway:** Malawi 1997 and Belize 1993/1994/1997 remain, by a wide
margin, the most robust findings — they top both rankings and are flagged by
both methods with a very large gap over the rest of the list. Under the
corrected (S-estimator, unweighted) methodology, Mozambique 2014, Burkina
Faso 1998, Mozambique 2002, Central African Republic 2008 and Rwanda 2000
drop out of the "always problematic" list (they remain LIS-flagged but no
longer clear the tighter, less-sensitive S-estimator threshold), while Rwanda
2013 (EICV-IV) newly joins it. The list is now smaller (10 vs. the earlier
14) and slightly more concentrated in Sub-Saharan African consumption
surveys, with Seychelles 2013 (income) as the one exception, reinforcing
that this is not purely a consumption-vs-income artifact.

See `outputs/always_problematic_surveys.csv` for the full list with values,
and `outputs/lis_vs_outdetect_comparison.csv` for the full 1,825-row
side-by-side comparison (includes `welfare_type`, consumption vs. income,
merged from `interpolated_means.csv`).
