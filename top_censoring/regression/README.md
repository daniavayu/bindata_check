# Regression analyses

Reserved for the Malawi 1997 regression analysis using household and demographic characteristics.

The planned comparison is between observed welfare and welfare predicted from observable household characteristics. This analysis is kept separate from the distributional analyses in `pareto/` and `gb2/`.

## Malawi 1997 implementation

The Malawi 1997 country file does not contain roof material, water access, or
education variables. The implemented model therefore predicts log welfare in
2021 PPP USD per day from the available variables: urban status, log household
size, and subnational area (`subnatid`). The spatial/temporal food and non-food
indices are not used as predictors because they are price-adjustment fields,
not household characteristics.

The model uses the survey weights in a weighted least-squares fit. Predictions
are generated with deterministic five-fold cross-fitting so that each
household is compared with a prediction from a model that did not train on
that household. Households with an absolute standardized cross-fitted log
residual of at least 3 are flagged for review. A flag means that reported
welfare is unusually high or low relative to the included characteristics; it
does not establish that the report is erroneous.

Run from the `bindata_check` directory:

```powershell
python top_censoring\regression\01_regression_mwi_1997.py
```

Outputs are written to `top_censoring/regression/outputs/`, including the
household predictions and flags, model summary, fitted coefficients, and an
observed-versus-predicted plot.
