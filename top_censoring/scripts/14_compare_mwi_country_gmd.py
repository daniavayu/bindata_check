from pathlib import Path

import numpy as np
import pandas as pd


root = Path(__file__).resolve().parents[2]
country_path = root / "01-input" / "country" / "MWI_1997.dta"
gmd_path = root / "01-input" / "GMD_all_2017.dta"

country = pd.read_stata(country_path)
country["welfare_2021"] = country["welfare"] / (
    country["cpi2021"] * country["icp2021"] * 365
)
country["welfare_2017"] = country["welfare"] / (
    country["cpi2017"] * country["icp2017"] * 365
)
country_valid = country.loc[
    country["welfare_2021"].gt(0)
    & np.isfinite(country["welfare_2021"])
    & country["weight"].gt(0)
    & np.isfinite(country["weight"])
].copy()

parts = []
reader = pd.read_stata(
    gmd_path,
    iterator=True,
    chunksize=100_000,
    columns=[
        "welfare_dec",
        "weight",
        "code",
        "year",
        "use_microdata",
        "use_bin",
        "use_groupdata",
    ],
)
for chunk in reader:
    mask = chunk["code"].astype(str).str.upper().eq("MWI")
    mask = mask & pd.to_numeric(chunk["year"], errors="coerce").eq(1997)
    mask = mask & pd.to_numeric(chunk["use_microdata"], errors="coerce").eq(1)
    mask = mask & pd.to_numeric(chunk["use_bin"], errors="coerce").eq(0)
    mask = mask & pd.to_numeric(chunk["use_groupdata"], errors="coerce").eq(0)
    if mask.any():
        parts.append(chunk.loc[mask, ["welfare_dec", "weight"]])

gmd = pd.concat(parts, ignore_index=True)
gmd["welfare_dec"] = pd.to_numeric(gmd["welfare_dec"], errors="coerce")
gmd["weight"] = pd.to_numeric(gmd["weight"], errors="coerce")
gmd_valid = gmd.loc[
    gmd["welfare_dec"].gt(0)
    & np.isfinite(gmd["welfare_dec"])
    & gmd["weight"].gt(0)
    & np.isfinite(gmd["weight"])
].copy()

country_values = country_valid["welfare_2021"].to_numpy(dtype=float)
country_values_2017 = country_valid["welfare_2017"].to_numpy(dtype=float)
country_weights = country_valid["weight"].to_numpy(dtype=float)
gmd_values = gmd_valid["welfare_dec"].to_numpy(dtype=float)
gmd_weights = gmd_valid["weight"].to_numpy(dtype=float)

ppp_factor = np.median(country_values / country_values_2017)
country_sorted = np.sort(country_values)
gmd_scaled_sorted = np.sort(gmd_values * ppp_factor)
common = min(len(country_sorted), len(gmd_scaled_sorted))
welfare_differences = np.abs(country_sorted[:common] - gmd_scaled_sorted[:common])

country_weights_sorted = np.sort(country_weights)
gmd_weights_sorted = np.sort(gmd_weights)
weight_differences = np.abs(country_weights_sorted[:common] - gmd_weights_sorted[:common])

country_mean = np.sum(country_values * country_weights) / country_weights.sum()
gmd_mean_2017 = np.sum(gmd_values * gmd_weights) / gmd_weights.sum()
gmd_mean_2021 = gmd_mean_2017 * ppp_factor

print(f"Country valid observations: {len(country_values):,}")
print(f"GMD valid observations: {len(gmd_values):,}")
print(f"PPP 2017-to-2021 factor: {ppp_factor:.12f}")
print(f"Welfare max absolute difference after scaling: {welfare_differences.max():.12f}")
print(f"Welfare mean absolute difference after scaling: {welfare_differences.mean():.12f}")
print(f"Welfare values equal within 1e-8: {np.allclose(country_sorted[:common], gmd_scaled_sorted[:common], rtol=1e-8, atol=1e-8)}")
print(f"Country weight sum: {country_weights.sum():.12f}")
print(f"GMD weight sum: {gmd_weights.sum():.12f}")
print(f"Weight max absolute difference: {weight_differences.max():.12f}")
print(f"Weight mean absolute difference: {weight_differences.mean():.12f}")
print(f"Weights equal within 1e-8: {np.allclose(country_weights_sorted[:common], gmd_weights_sorted[:common], rtol=1e-8, atol=1e-8)}")
print(f"Country weighted mean (PPP 2021): {country_mean:.12f}")
print(f"GMD weighted mean after PPP conversion: {gmd_mean_2021:.12f}")