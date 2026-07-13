from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/Users/danielaavayu/Desktop/World Bank/Bottom Coding")
INPUT = ROOT / "01-input"
YUKI = INPUT / "nuevosyuki"

fillgaps_path = YUKI / "fillgaps.dta"
bins_path = YUKI / "GlobalDist1000bins_1990_2026_20260324_2021_01_02_PROD.dta"
out_path = YUKI / "pip_1kbins_fillgaps.dta"

fillgaps = pd.read_stata(fillgaps_path, convert_categoricals=False)

chunks = []

for chunk in pd.read_stata(
    bins_path,
    columns=["code", "year", "quantile", "welf", "pop"],
    chunksize=500_000,
):
    chunk["hc_30"] = chunk["welf"] < 3

    collapsed = (
        chunk.groupby(["code", "year"], as_index=False)
        .apply(lambda g: pd.Series({
            "mean_1kbins": np.average(g["welf"], weights=g["pop"]),
            "headcount_1kbins": np.average(g["hc_30"], weights=g["pop"]),
            "pop": g["pop"].sum(),
        }))
        .reset_index(drop=True)
    )

    chunks.append(collapsed)

bins_collapsed = (
    pd.concat(chunks, ignore_index=True)
    .groupby(["code", "year"], as_index=False)
    .apply(lambda g: pd.Series({
        "mean_1kbins": np.average(g["mean_1kbins"], weights=g["pop"]),
        "headcount_1kbins": np.average(g["headcount_1kbins"], weights=g["pop"]),
        "pop": g["pop"].sum(),
    }))
    .reset_index(drop=True)
)

bins_collapsed = bins_collapsed.rename(columns={"code": "country_code"})

merged = fillgaps.merge(
    bins_collapsed,
    on=["country_code", "year"],
    how="inner",
)

merged.to_stata(out_path, write_index=False)

merged.head()