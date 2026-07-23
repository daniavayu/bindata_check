# Data Setup (Local Only)

This project expects local data folders that are not tracked in git.

## Expected Local Paths

- `01-input/`
- `02-output/`
- `work/`
- `interpolation_excercise/`

## Core Input Files Used by Notebooks

- `01-input/pip_1kbins.dta`
- `01-input/GlobalDist1000bins_1990_2026_20260324_2021_01_02_PROD.dta`
- `01-input/country/MWI_1997.dta`
- `01-input/country/MWI_2004.dta`
- `01-input/20260922/lineup/fillgaps.dta`
- `01-input/20260922/lineup/GlobalDist1000bins_1990_2026_20260922_2021_01_02_PROD.dta`
- `01-input/interpolated_means.dta`

## Why Data Is Excluded

GitHub rejects files over 100 MB, and several project inputs are much larger.
The repository tracks code and documentation; data is maintained locally.