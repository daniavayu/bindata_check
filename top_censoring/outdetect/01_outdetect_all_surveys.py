"""outdetect-style top-tail outlier screening for every GMD survey.

Replicates (in Python, since Stata's `outdetect` is not available here) the
robust z-score outlier identifier requested by the team:

    outdetect welfare, norm(ln) zscore(median, s) alpha(3) out(top)

- norm(ln)        -> work in log(welfare) space.
- zscore(median,s) -> center = weighted median of log(welfare);
                      scale  = weighted MAD of log(welfare), scaled by 1.4826
                      so it is comparable to a standard deviation under
                      normality (this is the classical Hampel-identifier
                      scale). NOTE: the exact letter "s" in the original
                      Stata syntax could not be confirmed against the help
                      file (no internet access in this environment); MAD is
                      documented here as the working assumption. See README.md.
- alpha(3)        -> flag values more than 3 robust "scales" above the median.
- out(top)        -> only the upper tail is flagged (top-tail outliers).

Run on the same universe as the existing LIS-ceiling screening
(top_censoring/scripts/15_gmd_all_surveys_lis_indicator_comparison.py):
every (country, year, survey) group with usable microdata in GMD_all_2017.dta.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUTPUT_PATH = Path(__file__).resolve().parent / "outputs" / "outdetect_all_surveys_summary.csv"

ALPHA = 3
MAD_SCALE_FACTOR = 1.4826  # makes MAD comparable to SD under normality


def weighted_quantile(values, weights, probability):
    order = np.argsort(values)
    sorted_values = np.asarray(values)[order]
    sorted_weights = np.asarray(weights)[order]
    cumulative = np.cumsum(sorted_weights) / sorted_weights.sum()
    return sorted_values[np.searchsorted(cumulative, probability, side="left")]


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

for (survey, country, year), data in microdata.groupby(["Survey", "Country", "Year"], sort=True):
    values = data["welfare"].to_numpy(dtype=float)
    weights = data["weight"].to_numpy(dtype=float)
    log_values = np.log(values)

    median_log = weighted_quantile(log_values, weights, 0.5)
    mad_log = weighted_quantile(np.abs(log_values - median_log), weights, 0.5)
    scale_log = MAD_SCALE_FACTOR * mad_log
    threshold = np.exp(median_log + ALPHA * scale_log) if scale_log > 0 else np.inf

    top_outlier = values > threshold
    capped_values = np.minimum(values, threshold)
    original = indicators(values, weights)
    capped = indicators(capped_values, weights)

    result = {
        "Country": country,
        "Year": year,
        "Survey": survey,
        "Valid records": len(values),
        "Median log welfare": median_log,
        "MAD log welfare": mad_log,
        "Robust scale (1.4826*MAD)": scale_log,
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
