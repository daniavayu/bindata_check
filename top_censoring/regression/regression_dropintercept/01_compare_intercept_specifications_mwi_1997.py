"""Compare weighted cross-fitted welfare regressions with and without an intercept.

Both models predict log welfare from urban status, log household size, and
subnational area. The no-intercept model includes all categories of the
categorical predictors, so no reference category is silently forced to zero.
"""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
INPUT_PATH = ROOT / "01-input" / "country" / "MWI_1997.dta"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
N_FOLDS = 5
FLAG_THRESHOLD = 3.0
TOP_N = 5


def load_data():
    columns = [
        "hhid",
        "welfare",
        "weight",
        "cpi2021",
        "icp2021",
        "urban",
        "hsize",
        "subnatid",
    ]
    data = pd.read_stata(INPUT_PATH, columns=columns)
    data["welfare_ppp_day"] = data["welfare"] / (
        data["cpi2021"] * data["icp2021"] * 365
    )
    valid = (
        np.isfinite(data["welfare_ppp_day"])
        & data["welfare_ppp_day"].gt(0)
        & np.isfinite(data["weight"])
        & data["weight"].gt(0)
        & np.isfinite(data["hsize"])
        & data["urban"].notna()
        & data["subnatid"].notna()
    )
    return data.loc[valid].reset_index(drop=True)


def make_design_matrices(data):
    predictors = data[["urban", "hsize", "subnatid"]].copy()
    predictors["log_hsize"] = np.log1p(predictors.pop("hsize"))
    categorical = pd.get_dummies(
        predictors[["urban", "subnatid"]],
        columns=["urban", "subnatid"],
        drop_first=False,
        dtype=float,
    )
    no_intercept = categorical.join(predictors[["log_hsize"]])
    with_intercept = no_intercept.drop(
        columns=[column for column in no_intercept.columns if column.startswith(("urban_", "subnatid_"))]
    )
    # Rebuild the reference-category specification used by the original script.
    reference = pd.get_dummies(
        predictors[["urban", "subnatid"]],
        columns=["urban", "subnatid"],
        drop_first=True,
        dtype=float,
    ).join(predictors[["log_hsize"]])
    return (
        np.column_stack([np.ones(len(reference)), reference.to_numpy(float)]),
        ["intercept", *reference.columns.tolist()],
        no_intercept.to_numpy(float),
        no_intercept.columns.tolist(),
    )


def weighted_fit(design, target, weights):
    sqrt_weights = np.sqrt(weights / np.mean(weights))
    coefficients, _, _, _ = np.linalg.lstsq(
        design * sqrt_weights[:, None], target * sqrt_weights, rcond=None
    )
    return coefficients


def weighted_mean(values, weights):
    return np.sum(values * weights) / np.sum(weights)


def weighted_sd(values, weights):
    mean = weighted_mean(values, weights)
    return np.sqrt(weighted_mean((values - mean) ** 2, weights))


def weighted_r2(target, prediction, weights):
    centered = target - weighted_mean(target, weights)
    return 1 - np.sum(weights * (target - prediction) ** 2) / np.sum(weights * centered**2)


def cross_fitted_predictions(design, target, weights):
    predictions = np.full(len(target), np.nan)
    folds = np.arange(len(target)) % N_FOLDS
    for fold in range(N_FOLDS):
        train = folds != fold
        test = ~train
        coefficients = weighted_fit(design[train], target[train], weights[train])
        predictions[test] = design[test] @ coefficients
    return predictions


def model_results(design, target, weights):
    crossfit_log = cross_fitted_predictions(design, target, weights)
    residual_log = target - crossfit_log
    residual_scale = weighted_sd(residual_log, weights)
    return {
        "crossfit_log": crossfit_log,
        "predicted_welfare": np.exp(crossfit_log),
        "residual_log": residual_log,
        "standardized_residual": residual_log / residual_scale,
        "r2": weighted_r2(target, crossfit_log, weights),
        "rmse": np.sqrt(weighted_mean(residual_log**2, weights)),
        "residual_scale": residual_scale,
    }


def main():
    data = load_data()
    target = np.log(data["welfare_ppp_day"].to_numpy(float))
    weights = data["weight"].to_numpy(float)
    design, names, no_intercept_design, no_intercept_names = make_design_matrices(data)

    with_intercept = model_results(design, target, weights)
    without_intercept = model_results(no_intercept_design, target, weights)

    data["fitted_welfare_with_intercept"] = with_intercept["predicted_welfare"]
    data["fitted_welfare_without_intercept"] = without_intercept["predicted_welfare"]
    data["log_observed_welfare"] = np.log(data["welfare_ppp_day"])
    data["log_fitted_welfare"] = np.log(data["fitted_welfare_with_intercept"])
    data["ratio_with_intercept"] = data["welfare_ppp_day"] / data["fitted_welfare_with_intercept"]
    data["ratio_without_intercept"] = data["welfare_ppp_day"] / data["fitted_welfare_without_intercept"]
    data["flag_with_intercept"] = with_intercept["standardized_residual"] >= FLAG_THRESHOLD
    data["flag_without_intercept"] = without_intercept["standardized_residual"] >= FLAG_THRESHOLD

    top_columns = [
        "hhid",
        "welfare_ppp_day",
        "fitted_welfare_with_intercept",
        "fitted_welfare_without_intercept",
        "log_observed_welfare",
        "log_fitted_welfare",
        "ratio_with_intercept",
        "ratio_without_intercept",
        "urban",
        "hsize",
        "subnatid",
        "flag_with_intercept",
        "flag_without_intercept",
    ]
    top_table = data.nlargest(TOP_N, "welfare_ppp_day")[top_columns]

    summary = pd.DataFrame(
        [
            {
                "model": "with_intercept",
                "records": len(data),
                "predictors": "intercept + urban + log1p(hsize) + subnatid (reference categories)",
                "weighted_r2_crossfit_log_welfare": with_intercept["r2"],
                "weighted_rmse_crossfit_log_welfare": with_intercept["rmse"],
                "weighted_residual_sd_crossfit_log_welfare": with_intercept["residual_scale"],
                "flagged_high_records": int((with_intercept["standardized_residual"] >= FLAG_THRESHOLD).sum()),
            },
            {
                "model": "without_intercept",
                "records": len(data),
                "predictors": "all urban + all subnatid categories + log1p(hsize)",
                "weighted_r2_crossfit_log_welfare": without_intercept["r2"],
                "weighted_rmse_crossfit_log_welfare": without_intercept["rmse"],
                "weighted_residual_sd_crossfit_log_welfare": without_intercept["residual_scale"],
                "flagged_high_records": int((without_intercept["standardized_residual"] >= FLAG_THRESHOLD).sum()),
            },
        ]
    )

    coefficients = pd.concat(
        [
            pd.DataFrame({"model": "with_intercept", "term": names, "coefficient_log_welfare": weighted_fit(design, target, weights)}),
            pd.DataFrame({"model": "without_intercept", "term": no_intercept_names, "coefficient_log_welfare": weighted_fit(no_intercept_design, target, weights)}),
        ],
        ignore_index=True,
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    top_table.to_csv(OUTPUT_DIR / "top_5_observations_both_models.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "regression_model_comparison.csv", index=False)
    coefficients.to_csv(OUTPUT_DIR / "regression_coefficients_both_models.csv", index=False)
    data.to_csv(OUTPUT_DIR / "regression_predictions_both_models.csv", index=False)

    print(summary.to_string(index=False))
    print(f"\nTop {TOP_N} observations by observed welfare:")
    print(top_table.to_string(index=False))
    print(f"\nWrote comparison outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
