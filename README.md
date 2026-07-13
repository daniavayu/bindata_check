# Bottom Censoring

Replication and diagnostic analysis of World Bank 1,000-bin welfare distributions, focused on Malawi, to trace mean gaps between survey microdata, fillgaps interpolation, and GlobalDist top-tail adjustments.

## Repository Layout

- `00-code/`: legacy production scripts (`.do` and python helper).
- `docs/literature/`: reference papers and background documents.
- Root notebooks (`*.ipynb`): analysis and replication notebooks. These stay at root because many use `Path(".")` and expect data folders relative to the project root.
- Root scripts (`*.R`, `main.py`): utility and prototype scripts.

## Data Policy

Raw and heavy data are intentionally excluded from git:

- `01-input/`
- `02-output/`
- `work/`
- `interpolation_excercise/`
- large binary formats (`.dta`, `.fst`, `.parquet`, etc.)

See `docs/DATA_SETUP.md` for local data placement.

## Quick Start

1. Create and activate your Python environment.
2. Install dependencies from `pyproject.toml`.
3. Place required local data under `01-input/` (not tracked by git).
4. Open and run notebooks from project root in VS Code/Jupyter.

## Notes

- Branch `master` is configured with a cleaned history (no large raw data files).
- If you add new large datasets locally, `.gitignore` will keep them out of commits.
