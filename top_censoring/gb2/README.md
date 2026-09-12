# GB2 analyses

This folder contains the Malawi 1997 GB2 distribution-fit analysis.

The planned comparison is:

- GB2 fitted to the full positive welfare sample.
- GB2 fitted after removing the largest 1, 2, 5, 10, and 20 observations.
- Observed versus fitted survival probabilities at the 30 largest observed values.
- A log-log survival plot comparing all fitted GB2 curves with the empirical tail.

Run from the `bindata_check` directory:

```powershell
C:\WBG\Python313\python.exe top_censoring\gb2\01_fit_gb2_mwi_1997.py
```

The fit uses weighted maximum likelihood and the same Malawi country file and
2021 PPP conversion used by the existing Pareto analysis. Outputs are written
to `top_censoring/gb2/outputs/`.

## Malawi 1997 results

All six fits converged. The full-sample estimates are `a=4.100`,
`b=1.705`, `p=0.862`, and `q=0.407`. Removing the top 20 records changes
these to `a=3.855`, `b=1.707`, `p=0.930`, and `q=0.447`, indicating sensitivity
of the fitted upper tail to the extreme observations.

For the full-sample fit, the second-largest observed value has an empirical
survival probability about 166 times larger than the fitted GB2 survival
probability. The largest observation is even more extreme; its fitted
survival probability underflows to zero in double precision. These results
support treating the very top observations as influential relative to a smooth
GB2 tail, but they do not establish that the observations are invalid.

The model is diagnostic: excluding extreme observations and refitting the GB2
shows how much the estimated distribution changes, but does not by itself
prove that any observation is erroneous.
