from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUTPUT_PATH = ROOT / "top_censoring" / "outputs" / "all_surveys_four_indicators_lis_comparison.csv"


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
    q1 = weighted_quantile(np.log(values), weights, 0.25)
    q3 = weighted_quantile(np.log(values), weights, 0.75)
    lis_ceiling = np.exp(q3 + 3 * (q3 - q1))
    capped_values = np.minimum(values, lis_ceiling)
    original = indicators(values, weights)
    capped = indicators(capped_values, weights)

    result = {
        "Country": country,
        "Year": year,
        "Survey": survey,
        "Valid records": len(values),
        "LIS ceiling (PPP/day)": lis_ceiling,
        "Welfare share above ceiling (%)": 100
        * np.sum(values[values > lis_ceiling] * weights[values > lis_ceiling])
        / np.sum(values * weights),
    }
    for name in original:
        result[f"{name} original"] = original[name]
        result[f"{name} LIS-capped"] = capped[name]
        result[f"{name} change (%)"] = 100 * (capped[name] / original[name] - 1)
    results.append(result)

comparison = pd.DataFrame(results)
comparison = comparison.sort_values(["Country", "Year", "Survey"]).reset_index(drop=True)
comparison.to_csv(OUTPUT_PATH, index=False)
print(f"Wrote {len(comparison)} survey comparisons to {OUTPUT_PATH}")
print(comparison.head(10).to_string(index=False))