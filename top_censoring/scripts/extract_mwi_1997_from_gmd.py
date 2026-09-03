from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GMD_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUT_PATH = ROOT / "top_censoring" / "outputs" / "MWI_1997_from_GMD.csv"

parts = []
reader = pd.read_stata(GMD_PATH, iterator=True, chunksize=100_000)
for chunk in reader:
    year = pd.to_numeric(chunk["year"], errors="coerce")
    mask = chunk["code"].astype(str).str.upper().eq("MWI") & year.eq(1997)
    if mask.any():
        parts.append(chunk.loc[mask].copy())

if parts:
    output = pd.concat(parts, ignore_index=True)
else:
    output = pd.DataFrame()

output.to_csv(OUT_PATH, index=False)

print(f"rows={len(output)}")
print(f"cols={len(output.columns)}")
print(f"saved={OUT_PATH}")
