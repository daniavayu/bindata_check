"""Executive summary: before/after comparison for the surveys flagged by both methods.

Recomputes, directly from GMD_all_2017.dta, both the LIS (Tukey/IQR) ceiling
and the outdetect (median + S-estimator) threshold for the 10 survey-years
flagged as problematic by BOTH methods (see ../README.md, section 6), plus a
record-level deep dive on the flagship case, Malawi 1997 (IHS-I).

Outputs (clean, wide tables ready to paste into a slide or open in Excel):
    outputs/executive_summary_always_problematic.csv
  outputs/executive_summary_mwi_1997_flagship_top_values.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
ALWAYS_PROBLEMATIC_PATH = OUT_DIR / "always_problematic_surveys.csv"

ALPHA = 3
S_SCALE_FACTOR = 1.1926
S_ESTIMATOR_SUBSAMPLE_SIZE = 2000
RANDOM_SEED = 42

ALWAYS_PROBLEMATIC = list(
    pd.read_csv(ALWAYS_PROBLEMATIC_PATH)[["Country", "Year", "Survey"]]
    .itertuples(index=False, name=None)
)
FLAGSHIP = ("MWI", 1997, "IHS-I")


def s_estimator(values, rng):
    if len(values) > S_ESTIMATOR_SUBSAMPLE_SIZE:
        values = rng.choice(values, size=S_ESTIMATOR_SUBSAMPLE_SIZE, replace=False)
    diffs = np.abs(values[:, None] - values[None, :])
    row_medians = np.median(diffs, axis=1)
    return S_SCALE_FACTOR * np.median(row_medians)


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
        "Mean (PPP/day)": mean_welfare,
        "Gini": weighted_gini(values, weights),
        "MLD": mld,
        "Atkinson(2)": atkinson_2,
    }


# --- Load only the microdata for the 10 target surveys -----------------------------------
wanted = {(c, y, s) for c, y, s in ALWAYS_PROBLEMATIC}
columns = ["welfare_dec", "weight", "code", "year", "survname", "use_microdata", "use_bin", "use_groupdata"]
parts = []
reader = pd.read_stata(DATA_PATH, iterator=True, chunksize=100_000, columns=columns)
for chunk in reader:
    year = pd.to_numeric(chunk["year"], errors="coerce")
    country = chunk["code"].astype(str).str.upper()
    survey = chunk["survname"].astype(str)
    welfare = pd.to_numeric(chunk["welfare_dec"], errors="coerce")
    weight = pd.to_numeric(chunk["weight"], errors="coerce")
    microdata_flag = pd.to_numeric(chunk["use_microdata"], errors="coerce").eq(1)
    not_grouped = (
        pd.to_numeric(chunk["use_bin"], errors="coerce").eq(0)
        & pd.to_numeric(chunk["use_groupdata"], errors="coerce").eq(0)
    )
    valid = (
        microdata_flag
        & not_grouped
        & welfare.gt(0)
        & np.isfinite(welfare)
        & weight.gt(0)
        & np.isfinite(weight)
    )
    if valid.any():
        keep = pd.Series(
            list(zip(country.loc[valid], year.loc[valid].astype(int), survey.loc[valid])),
            index=welfare.loc[valid].index,
        ).isin(wanted)
        idx = welfare.loc[valid].index[keep]
        if len(idx):
            parts.append(
                pd.DataFrame(
                    {
                        "Country": country.loc[idx].to_numpy(),
                        "Year": year.loc[idx].astype(int).to_numpy(),
                        "Survey": survey.loc[idx].to_numpy(),
                        "welfare": welfare.loc[idx].to_numpy(),
                        "weight": weight.loc[idx].to_numpy(),
                    }
                )
            )

microdata = pd.concat(parts, ignore_index=True)
rng = np.random.default_rng(RANDOM_SEED)

# --- Table 1: boss-ready wide comparison for the 10 surveys -------------------------------
rows = []
flagship_values = None
flagship_weights = None
flagship_thresholds = None

for country, year, survey in ALWAYS_PROBLEMATIC:
    data = microdata[
        (microdata["Country"] == country) & (microdata["Year"] == year) & (microdata["Survey"] == survey)
    ]
    values = data["welfare"].to_numpy(dtype=float)
    weights = data["weight"].to_numpy(dtype=float)
    log_values = np.log(values)

    # LIS (Tukey/IQR) ceiling.
    q1 = weighted_quantile(log_values, weights, 0.25)
    q3 = weighted_quantile(log_values, weights, 0.75)
    lis_ceiling = np.exp(q3 + 3 * (q3 - q1))
    lis_flag = values > lis_ceiling
    lis_capped = np.minimum(values, lis_ceiling)

    # outdetect (median + S-estimator) threshold.
    median_log = np.median(log_values)
    s_log = s_estimator(log_values, rng)
    outdetect_threshold = np.exp(median_log + ALPHA * s_log)
    outdetect_flag = values > outdetect_threshold
    outdetect_capped = np.minimum(values, outdetect_threshold)

    raw = indicators(values, weights)
    lis_ind = indicators(lis_capped, weights)
    outdetect_ind = indicators(outdetect_capped, weights)

    row = {
        "Country": country,
        "Year": year,
        "Survey": survey,
        "Valid records": len(values),
        "LIS ceiling (PPP/day)": lis_ceiling,
        "outdetect threshold (PPP/day)": outdetect_threshold,
        "Population affected (%) - LIS": 100 * np.sum(weights[lis_flag]) / np.sum(weights),
        "Population affected (%) - outdetect": 100 * np.sum(weights[outdetect_flag]) / np.sum(weights),
        "Welfare affected (%) - LIS": 100 * np.sum(values[lis_flag] * weights[lis_flag]) / np.sum(values * weights),
        "Welfare affected (%) - outdetect": 100
        * np.sum(values[outdetect_flag] * weights[outdetect_flag])
        / np.sum(values * weights),
    }
    for name in raw:
        row[f"{name} - Raw"] = raw[name]
        row[f"{name} - LIS-capped"] = lis_ind[name]
        row[f"{name} - Change (%) LIS"] = 100 * (lis_ind[name] / raw[name] - 1)
        row[f"{name} - outdetect-capped"] = outdetect_ind[name]
        row[f"{name} - Change (%) outdetect"] = 100 * (outdetect_ind[name] / raw[name] - 1)
    rows.append(row)

    if (country, year, survey) == FLAGSHIP:
        flagship_values = values
        flagship_weights = weights
        flagship_thresholds = (lis_ceiling, outdetect_threshold, lis_flag, outdetect_flag)

summary = pd.DataFrame(rows)
OUT_DIR.mkdir(parents=True, exist_ok=True)
summary_path = OUT_DIR / "executive_summary_always_problematic.csv"
summary.to_csv(summary_path, index=False)
print(f"Wrote {len(summary)} rows to {summary_path}")

# --- Table 2: Malawi 1997 flagship deep dive (top 30 largest welfare values) -------------
lis_ceiling, outdetect_threshold, lis_flag, outdetect_flag = flagship_thresholds
order = np.argsort(flagship_values)[::-1][:30]
flagship_table = pd.DataFrame(
    {
        "Rank (largest first)": np.arange(1, len(order) + 1),
        "Welfare (PPP/day)": flagship_values[order],
        "Sampling weight": flagship_weights[order],
        "Flagged by LIS (> ceiling)": lis_flag[order],
        "Flagged by outdetect (> threshold)": outdetect_flag[order],
        "Flagged by both": lis_flag[order] & outdetect_flag[order],
    }
)
flagship_table.attrs["lis_ceiling"] = lis_ceiling
flagship_table.attrs["outdetect_threshold"] = outdetect_threshold
flagship_path = OUT_DIR / "executive_summary_mwi_1997_flagship_top_values.csv"
with open(flagship_path, "w", newline="") as f:
    f.write(f"# LIS ceiling (PPP/day): {lis_ceiling:.4f}\n")
    f.write(f"# outdetect threshold (PPP/day): {outdetect_threshold:.4f}\n")
    flagship_table.to_csv(f, index=False)
print(f"Wrote flagship top-30 table to {flagship_path}")
print(summary[["Country", "Year", "Survey", "Population affected (%) - LIS", "Population affected (%) - outdetect"]].to_string(index=False))
