"""outdetect-style top-tail outlier screening for every GMD survey.

Replicates (in Python, since Stata's `outdetect` is not available here) the
robust z-score outlier identifier requested by the team:

    outdetect welfare, norm(ln) zscore(median, s) alpha(3) out(top)

Methodology confirmed against Belotti, Mancini and Vecchi, "Outlier Detection
for Welfare Analysis" (World Bank Policy Research Working Paper):

- norm(ln)         -> work in log(welfare) space.
- zscore(median,s) -> center = UNWEIGHTED median of log(welfare);
                      scale  = UNWEIGHTED S-estimator (Rousseeuw & Croux,
                      1993) of log(welfare): S = 1.1926 * med_i{ med_j |y_i -
                      y_j| }. This is a *different* estimator from the MAD
                      (Hampel, 1974) despite both being median-based; the
                      paper's Table 2 lists them as separate options ("mad"
                      vs. "s"), and the requested syntax `zscore(median, s)`
                      maps to the S row, not the MAD row.
                      The paper explicitly notes center/scale are computed
                      UNWEIGHTED ("the weighted median does not... qualify as
                      a robust statistic"), even though the final indicator
                      comparison (Raw vs. Trimmed/capped) does use survey
                      weights. Both conventions are followed here.
- alpha(3)         -> flag values more than 3 robust "scales" above the median.
- out(top)         -> only the upper tail is flagged (top-tail outliers).

Run on the same universe as the existing LIS-ceiling screening
(top_censoring/scripts/15_gmd_all_surveys_lis_indicator_comparison.py):
every (country, year, survey) group with usable microdata in GMD_all_2017.dta.

Computational note: the S-estimator is O(n^2) (median of per-observation
medians of pairwise absolute differences). Survey sizes here range up to ~5
million records, so an exact computation is infeasible for most surveys.
When a group has more than S_ESTIMATOR_SUBSAMPLE_SIZE records, the S-estimator
(only) is computed on a fixed-seed random subsample of that size; the median
center and the resulting threshold/capping/indicators are still computed on
the FULL survey. This is a standard practical compromise for robust scale
estimators on large samples; see README.md for details.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUTPUT_PATH = Path(__file__).resolve().parent / "outputs" / "outdetect_all_surveys_summary.csv"

ALPHA = 3
S_SCALE_FACTOR = 1.1926  # Rousseeuw & Croux (1993) consistency constant for S
S_ESTIMATOR_SUBSAMPLE_SIZE = 2000  # cap for the O(n^2) S-estimator computation
RANDOM_SEED = 42


def s_estimator(values, rng):
    """Unweighted S-estimator of scale (Rousseeuw & Croux, 1993): S = 1.1926 * med_i{med_j|y_i - y_j|}."""
    if len(values) > S_ESTIMATOR_SUBSAMPLE_SIZE:
        values = rng.choice(values, size=S_ESTIMATOR_SUBSAMPLE_SIZE, replace=False)
    diffs = np.abs(values[:, None] - values[None, :])
    row_medians = np.median(diffs, axis=1)
    return S_SCALE_FACTOR * np.median(row_medians)


def weighted_mean(values, weights):
    return np.sum(values * weights) / np.sum(weights)


def weighted_gini(values, weights):
    order = np.argsort(values)
    sorted_values = np.asarray(values)[order]
    sorted_weights = np.asarray(weights)[order]
    population = np.r_[0.0, np.cumsum(sorted_weights) / sorted_weights.sum()]
    welfare = np.r_[0.0, np.cumsum(sorted_values * sorted_weights) / np.sum(sorted_values * sorted_weights)]
    return 1 - np.sum(np.diff(population) * (welfare[:-1] + welfare[1:]))


def indicators(values, weights):
    mean_welfare = weighted_mean(values, weights)
    mld = np.log(mean_welfare) - weighted_mean(np.log(values), weights)
    atkinson_2 = 1 - (1 / weighted_mean(1 / values, weights)) / mean_welfare
    return {
        "mean_ppp_day": mean_welfare,
        "gini": weighted_gini(values, weights),
        "mld": mld,
        "atkinson_2": atkinson_2,
    }


columns = [
    "welfare_dec",
    "weight",
    "code",
    "year",
    "survname",
    "use_microdata",
    "use_bin",
    "use_groupdata",
]
parts = []

reader = pd.read_stata(DATA_PATH, iterator=True, chunksize=100_000, columns=columns)
for chunk in reader:
    year = pd.to_numeric(chunk["year"], errors="coerce")
    country = chunk["code"].astype(str).str.upper()
    survey = chunk["survname"].astype(str)

    welfare = pd.to_numeric(chunk["welfare_dec"], errors="coerce")
    weight = pd.to_numeric(chunk["weight"], errors="coerce")
    microdata = pd.to_numeric(chunk["use_microdata"], errors="coerce").eq(1)
    not_grouped = (
        pd.to_numeric(chunk["use_bin"], errors="coerce").eq(0)
        & pd.to_numeric(chunk["use_groupdata"], errors="coerce").eq(0)
    )
    valid = (
        microdata
        & not_grouped
        & welfare.gt(0)
        & np.isfinite(welfare)
        & weight.gt(0)
        & np.isfinite(weight)
    )
    if valid.any():
        parts.append(
            pd.DataFrame(
                {
                    "Survey": survey.loc[valid].to_numpy(),
                    "Country": country.loc[valid].to_numpy(),
                    "Year": year.loc[valid].astype(int).to_numpy(),
                    "welfare": welfare.loc[valid].to_numpy(),
                    "weight": weight.loc[valid].to_numpy(),
                }
            )
        )

microdata = pd.concat(parts, ignore_index=True)
results = []
rng = np.random.default_rng(RANDOM_SEED)

for (survey, country, year), data in microdata.groupby(["Survey", "Country", "Year"], sort=True):
    values = data["welfare"].to_numpy(dtype=float)
    weights = data["weight"].to_numpy(dtype=float)
    log_values = np.log(values)

    median_log = np.median(log_values)  # unweighted, per Belotti/Mancini/Vecchi
    s_log = s_estimator(log_values, rng)  # unweighted S-estimator (Rousseeuw & Croux, 1993)
    threshold = np.exp(median_log + ALPHA * s_log) if s_log > 0 else np.inf

    top_outlier = values > threshold
    capped_values = np.minimum(values, threshold)
    original = indicators(values, weights)
    capped = indicators(capped_values, weights)

    result = {
        "Country": country,
        "Year": year,
        "Survey": survey,
        "Valid records": len(values),
        "Median log welfare (unweighted)": median_log,
        "S estimate, log space (unweighted, 1.1926 scaled)": s_log,
        "S-estimator subsampled (Y/N)": "Y" if len(values) > S_ESTIMATOR_SUBSAMPLE_SIZE else "N",
        "outdetect threshold (PPP/day)": threshold,
        "Records flagged (top)": int(np.sum(top_outlier)),
        "Population share flagged (%)": 100 * np.sum(weights[top_outlier]) / np.sum(weights),
        "Welfare share flagged (%)": 100
        * np.sum(values[top_outlier] * weights[top_outlier])
        / np.sum(values * weights)
        if top_outlier.any()
        else 0.0,
    }
    for name in original:
        result[f"{name} original"] = original[name]
        result[f"{name} outdetect-capped"] = capped[name]
        result[f"{name} change (%)"] = 100 * (capped[name] / original[name] - 1)
    results.append(result)

summary = pd.DataFrame(results)
summary = summary.sort_values("Welfare share flagged (%)", ascending=False).reset_index(drop=True)
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(OUTPUT_PATH, index=False)
print(f"Wrote {len(summary)} survey rows to {OUTPUT_PATH}")
print(summary[["Country", "Year", "Survey", "Records flagged (top)", "Welfare share flagged (%)"]].head(30).to_string(index=False))
