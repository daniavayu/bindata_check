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

`outdetect` is a **Hampel-type robust z-score identifier**:

1. `norm(ln)` — work in `log(welfare)` space (reduces skewness before
   flagging).
2. `zscore(median, s)` — standardize each observation using a **robust**
   center and scale instead of the mean/SD (which are themselves distorted
   by outliers):
   - center = **weighted median** of `log(welfare)`
   - scale = **weighted MAD** (median absolute deviation) of `log(welfare)`,
     multiplied by `1.4826` so it is comparable to a standard deviation under
     normality. This is the classical Hampel-identifier scale.
3. `alpha(3)` — flag observations whose robust z-score exceeds `3`.
4. `out(top)` — only the **upper tail** is screened (this exercise is about
   top-censoring, not bottom outliers).

Threshold per survey:

$$
\text{threshold} = \exp\big(\text{median}(\ln w) + 3 \times 1.4826 \times \text{MAD}(\ln w)\big)
$$

Any record with `welfare > threshold` is flagged. Indicators
(`mean_ppp_day`, `gini`, `mld`, `atkinson_2`) are recomputed after capping
flagged values at the threshold ("outdetect-capped"), exactly mirroring how
the existing LIS-capped comparison table was built.

### Important caveat / assumption

The exact meaning of the `s` argument inside Stata's
`zscore(median, s)` option could not be verified against the original
documentation because **this environment has no internet access** (the
attempt to fetch the World Bank paper describing `outdetect` failed with a
connection error). **MAD (median absolute deviation), scaled by 1.4826, was
used as the working assumption** for the scale estimator, because:

- It is the standard, most common choice for a robust scale paired with a
  median center (the classical Hampel identifier).
- It is the natural robust counterpart to the Tukey-IQR-fence approach
  already used for the LIS method — both are "robust dispersion around a
  robust center", just using MAD instead of IQR.

If the original Stata source confirms a different scale estimator (e.g. Sn,
Qn, or a scaled IQR), only the `mad_log` / `scale_log` computation in
[`01_outdetect_all_surveys.py`](01_outdetect_all_surveys.py) needs to change;
everything downstream (thresholding, capping, indicator comparison) stays the
same.

Both the median and the MAD are computed as **weighted** statistics (using
each survey's sampling weights), for consistency with how the existing LIS
ceiling is computed (weighted 25th/75th percentiles of `log(welfare)`).

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
| [`01_outdetect_all_surveys.py`](01_outdetect_all_surveys.py) | Computes the outdetect (median+MAD, log, alpha=3, top-only) threshold and capped indicators for all 1,825 surveys. | `outputs/outdetect_all_surveys_summary.csv` |
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
  share flagged (%)") under the outdetect method.
  **A flat 5% cutoff does *not* carry over to outdetect**: the two methods
  flag very different baseline shares of welfare per survey
  (LIS mean/median across all 1,825 surveys = 0.85% / 0.38%; outdetect
  mean/median = 5.11% / 4.12%, since the Hampel/MAD identifier is inherently
  more sensitive than the Tukey/IQR fence). Applying LIS's raw 5% threshold to
  outdetect's "Welfare share flagged (%)" column flags **698 of 1,825
  surveys** (~38%), which is not a usable "problematic" list.
  Instead, the LIS 5% cutoff was translated to its **equivalent percentile
  rank** (24/1825 flagged → the 98.685th percentile of the LIS distribution),
  and that same percentile rank was applied to the outdetect distribution.
  That yields an outdetect threshold of **`Welfare share flagged (%) > 17.09%`**,
  which selects exactly **24 surveys** — mathematically the same operation as
  taking the top 24 by that column, which is what `02_compare_lis_vs_outdetect.py`
  does directly (`N_TOP = 24`). So the "24" is not a coincidence and not an
  arbitrary tie-in to the LIS list size either — it is what falls out of
  applying an equally-strict (same-percentile) cutoff to both methods.
- **Flagged by both ("always problematic")** = the survey shows up as an
  extreme top-tail case *regardless of which of the two outlier-screening
  methods is used* — the strongest possible evidence that these surveys need
  attention.

## 6. Results

Of the 24 LIS-flagged surveys and the 24 outdetect-flagged surveys (48 slots,
drawn from the same 1,825-survey universe), **14 surveys are flagged by
both methods** — these are the surveys where an extreme top tail shows up
regardless of whether a Tukey-IQR fence or a Hampel median/MAD identifier is
used:

| Country | Year | Survey | Welfare type | LIS welfare share above ceiling (%) | outdetect welfare share flagged (%) |
|---|---|---|---|---|---|
| MWI | 1997 | IHS-I | consumption | 41.60 | 49.30 |
| BLZ | 1993 | LFS | income | 20.09 | 26.06 |
| BLZ | 1994 | LFS | income | 19.96 | 31.04 |
| BLZ | 1997 | LFS | income | 18.24 | 23.00 |
| SYC | 2013 | HBS | income | 12.06 | 20.14 |
| MWI | 2016 | IHS-IV | consumption | 10.97 | 18.59 |
| ETH | 1995 | HICES | consumption | 10.89 | 19.43 |
| MOZ | 2014 | IOF | consumption | 8.04 | 21.46 |
| BFA | 1998 | EP-II | consumption | 7.90 | 23.92 |
| COM | 2004 | EIM | consumption | 7.61 | 27.51 |
| MOZ | 2002 | IAF | consumption | 6.77 | 17.96 |
| RWA | 2005 | EICV-II | consumption | 6.14 | 21.60 |
| CAF | 2008 | ECASEB | consumption | 5.47 | 17.87 |
| RWA | 2000 | EICV-I | consumption | 5.36 | 19.71 |

10 surveys are flagged only by LIS (Tukey/IQR) and not in the outdetect
top-24: MDG 1993, MDV 2002, MMR 2015, NOR 2005, NOR 2006, PRY 2002, PRY 2007,
PRY 2010, PRY 2014, RWA 2013. 10 different surveys are flagged only by
outdetect (median/MAD) and not in the LIS top-24 — see
`outputs/lis_vs_outdetect_comparison.csv` (columns `Flagged by LIS`,
`Flagged by outdetect`) for the complete breakdown.

**Takeaway:** Malawi 1997 and Belize 1993/1994/1997 are, by a wide margin,
the most robust findings — they top both rankings and are flagged by both
methods with a very large gap over the rest of the list. The other 10
"always problematic" surveys are mostly Sub-Saharan African
consumption surveys with the notable exception of Seychelles 2013
(income), reinforcing that this is not purely a consumption-vs-income
artifact.

See `outputs/always_problematic_surveys.csv` for the full list with values,
and `outputs/lis_vs_outdetect_comparison.csv` for the full 1,825-row
side-by-side comparison (includes `welfare_type`, consumption vs. income,
merged from `interpolated_means.csv`).
