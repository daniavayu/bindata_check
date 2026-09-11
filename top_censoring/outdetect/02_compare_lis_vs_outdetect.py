"""Compare the LIS Tukey-fence screening (existing method) against the
new outdetect-style (median + MAD, log space, alpha=3, top-tail only)
screening, over the same ~1825 survey-year universe from GMD_all_2017.dta.

Outputs:
  - lis_vs_outdetect_comparison.csv : full universe, both methods side by side,
    with flags for which method(s) singled out each survey as problematic.
  - always_problematic_surveys.csv : the subset flagged by BOTH methods
    (the direct answer to "which surveys are problematic no matter the
    outlier-detection method used").

Flagging rule (documented assumption, see README.md):
  - LIS-flagged  = survey is one of the 24 rows in
    24_surveys_four_indicators_lis_comparison.csv (the team's existing
    ground-truth "problematic" list).
  - outdetect-flagged = survey is among the top 24 surveys (by the same
    metric, "welfare share flagged/above ceiling (%)") under the outdetect
    method, i.e. the same list *size* as the LIS list, so the two methods
    are compared on equal footing rather than an arbitrarily different cutoff.
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]  # top_censoring/
OUT_DIR = Path(__file__).resolve().parent / "outputs"

LIS_FULL = BASE / "outputs" / "all_surveys_four_indicators_lis_comparison.csv"
LIS_TOP24 = BASE / "outputs" / "24_surveys_four_indicators_lis_comparison.csv"
OUTDETECT_FULL = OUT_DIR / "outdetect_all_surveys_summary.csv"
WELFARE_TYPE = (
    Path(__file__).resolve().parents[2]
    / "corrected_method"
    / "01-processing"
    / "raw"
    / "interpolated_means.csv"
)

N_TOP = 24  # matches the size of the existing ground-truth LIS list

lis_full = pd.read_csv(LIS_FULL)
lis_top24 = pd.read_csv(LIS_TOP24)
outdetect_full = pd.read_csv(OUTDETECT_FULL)

lis_top24_keys = set(zip(lis_top24["Country"], lis_top24["Year"]))

outdetect_top24 = (
    outdetect_full.sort_values("Welfare share flagged (%)", ascending=False)
    .head(N_TOP)
)
outdetect_top24_keys = set(zip(outdetect_top24["Country"], outdetect_top24["Year"]))

merged = lis_full.merge(
    outdetect_full,
    on=["Country", "Year", "Survey"],
    suffixes=(" (LIS)", " (outdetect)"),
    how="outer",
)

merged["Flagged by LIS"] = merged.apply(
    lambda r: (r["Country"], r["Year"]) in lis_top24_keys, axis=1
)
merged["Flagged by outdetect"] = merged.apply(
    lambda r: (r["Country"], r["Year"]) in outdetect_top24_keys, axis=1
)
merged["Flagged by both (always problematic)"] = (
    merged["Flagged by LIS"] & merged["Flagged by outdetect"]
)

welfare_type = pd.read_csv(WELFARE_TYPE, usecols=["country_code", "year", "welfare_type"])
welfare_type = welfare_type.drop_duplicates(subset=["country_code", "year"])
merged = merged.merge(
    welfare_type,
    left_on=["Country", "Year"],
    right_on=["country_code", "year"],
    how="left",
).drop(columns=["country_code", "year"])

merged = merged.sort_values(
    ["Flagged by both (always problematic)", "Welfare share above ceiling (%)"],
    ascending=[False, False],
)

OUT_DIR.mkdir(parents=True, exist_ok=True)
comparison_path = OUT_DIR / "lis_vs_outdetect_comparison.csv"
merged.to_csv(comparison_path, index=False)

always_problematic = merged[merged["Flagged by both (always problematic)"]]
always_problematic_path = OUT_DIR / "always_problematic_surveys.csv"
always_problematic.to_csv(always_problematic_path, index=False)

print(f"Wrote {len(merged)} rows to {comparison_path}")
print(f"Wrote {len(always_problematic)} rows to {always_problematic_path}")
print()
print("Surveys flagged as problematic by BOTH methods:")
print(
    always_problematic[
        ["Country", "Year", "Survey", "welfare_type", "Welfare share above ceiling (%)", "Welfare share flagged (%)"]
    ].to_string(index=False)
)
print()
print(f"LIS-only count: {(merged['Flagged by LIS'] & ~merged['Flagged by outdetect']).sum()}")
print(f"outdetect-only count: {(~merged['Flagged by LIS'] & merged['Flagged by outdetect']).sum()}")
