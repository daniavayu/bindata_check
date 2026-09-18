"""One-time fast export of the filtered microdata used by 01_outdetect_all_surveys.R.

R's haven::read_dta() has no chunked-reading mode, so pointing it directly at
the 3.6 GB GMD_all_2017.dta file (even with col_select) forces it to parse
the whole file's labelled-value structures in one pass, which is extremely
slow (tens of minutes) compared to pandas' chunked Stata reader.

This script reuses the exact same chunked read + filter logic as
../outdetect/01_outdetect_all_surveys.py to produce a small CSV containing
only the columns and rows the R script actually needs (valid microdata
records: use_microdata==1, use_bin==0, use_groupdata==0, welfare>0 finite,
weight>0 finite). 01_outdetect_all_surveys.R then reads this CSV (via
data.table::fread, seconds instead of tens of minutes) and does the actual
EXACT S-estimator computation (robustbase::Sn()) in R, per the user's
request. No statistics are computed here -- this is I/O only.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = ROOT / "01-input" / "GMD_all_2017.dta"
OUTPUT_PATH = Path(__file__).resolve().parent / "outputs" / "_microdata_filtered.csv"

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

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
if OUTPUT_PATH.exists():
    OUTPUT_PATH.unlink()

total_rows = 0
total_chunks = 0
reader = pd.read_stata(DATA_PATH, iterator=True, chunksize=100_000, columns=columns)
for chunk in reader:
    total_chunks += 1
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
        out = pd.DataFrame(
            {
                "Survey": survey.loc[valid].to_numpy(),
                "Country": country.loc[valid].to_numpy(),
                "Year": year.loc[valid].astype(int).to_numpy(),
                "welfare": welfare.loc[valid].to_numpy(),
                "weight": weight.loc[valid].to_numpy(),
            }
        )
        out.to_csv(OUTPUT_PATH, mode="a", header=(total_rows == 0), index=False)
        total_rows += len(out)

print(f"Read {total_chunks} chunks; wrote {total_rows:,} rows to {OUTPUT_PATH}")

