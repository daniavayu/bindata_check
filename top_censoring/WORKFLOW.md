# Malawi 1997 top-censoring workflow

Run the scripts from the repository root in this order:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" ".\top_censoring\scripts\01_lis_top_diagnostics.R"
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" ".\top_censoring\scripts\02_mwi_1997_pareto.R"
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" ".\top_censoring\scripts\03_mwi_1997_index_sensitivity.R"
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" ".\top_censoring\scripts\04_mwi_1997_dagum_sensitivity.R"
```

## What each script does

| Script | Question answered | Main output |
|---|---|---|
| `01_lis_top_diagnostics.R` | Which values exceed the LIS top ceiling? | LIS thresholds and shares affected |
| `02_mwi_1997_pareto.R` | Does a Pareto tail fitted from p90 describe the upper tail? | Pareto CCDF graph and parameter summary |
| `03_mwi_1997_index_sensitivity.R` | How influential are the two largest values for reported measures? | Mean, poverty, and inequality sensitivity table |
| `04_mwi_1997_dagum_sensitivity.R` | Does a full-distribution Dagum fit materially change after excluding them? | Side-by-side Dagum graph and fit summary |

## Interpretation order

1. LIS flags candidates; it does not prove an observation is wrong.
2. Pareto is fitted only to the upper tail and assesses whether very large
   observations are plausible under that tail.
3. The index table shows whether flagged observations matter in practice.
4. Dagum is fitted to the entire positive distribution. The fit without the two
   largest observations is a sensitivity exercise, not an automatic deletion
   rule.

All generated files are written to `top_censoring/outputs/`. No script modifies
the survey microdata in `01-input/country/MWI_1997.dta`.
