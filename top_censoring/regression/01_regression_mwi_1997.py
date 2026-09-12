"""Predict Malawi 1997 welfare from available household characteristics.

The Malawi 1997 country file does not contain roof, water, or education
variables. This analysis therefore uses the available household/geographic
predictors: urban status, household size, and subnational area.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = ROOT / "01-input" / "country" / "MWI_1997.dta"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
N_FOLDS = 5
FLAG_THRESHOLD = 3.0


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


def make_design_matrix(data):
    predictors = data[["urban", "hsize", "subnatid"]].copy()
    predictors["log_hsize"] = np.log1p(predictors.pop("hsize"))
    design = pd.get_dummies(
        predictors.drop(columns="log_hsize").join(predictors[["log_hsize"]]),
        columns=["urban", "subnatid"],
        drop_first=True,
        dtype=float,
    )
    return np.column_stack([np.ones(len(design)), design.to_numpy(float)]), [
        "intercept",
        *design.columns.tolist(),
    ]


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


def main():
    data = load_data()
    design, coefficient_names = make_design_matrix(data)
    target = np.log(data["welfare_ppp_day"].to_numpy(float))
    weights = data["weight"].to_numpy(float)

    full_coefficients = weighted_fit(design, target, weights)
    full_prediction_log = design @ full_coefficients

    crossfit_prediction_log = np.full(len(data), np.nan)
    folds = np.arange(len(data)) % N_FOLDS
    for fold in range(N_FOLDS):
        train = folds != fold
        test = ~train
        coefficients = weighted_fit(design[train], target[train], weights[train])
        crossfit_prediction_log[test] = design[test] @ coefficients

    residual_log = target - crossfit_prediction_log
    residual_scale = weighted_sd(residual_log, weights)
    data["predicted_welfare_ppp_day"] = np.exp(crossfit_prediction_log)
    data["in_sample_predicted_welfare_ppp_day"] = np.exp(full_prediction_log)
    data["log_residual_crossfit"] = residual_log
    data["actual_to_predicted_ratio"] = (
        data["welfare_ppp_day"] / data["predicted_welfare_ppp_day"]
    )
    data["standardized_log_residual"] = residual_log / residual_scale
    data["inconsistency_flag"] = (
        data["standardized_log_residual"].abs() >= FLAG_THRESHOLD
    )
    data["direction"] = np.where(
        data["standardized_log_residual"] >= FLAG_THRESHOLD,
        "reported welfare unusually high",
        np.where(
            data["standardized_log_residual"] <= -FLAG_THRESHOLD,
            "reported welfare unusually low",
            "",
        ),
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    prediction_columns = [
        "hhid",
        "subnatid",
        "urban",
        "hsize",
        "weight",
        "welfare_ppp_day",
        "predicted_welfare_ppp_day",
        "in_sample_predicted_welfare_ppp_day",
        "actual_to_predicted_ratio",
        "log_residual_crossfit",
        "standardized_log_residual",
        "inconsistency_flag",
        "direction",
    ]
    predictions = data[prediction_columns].sort_values(
        "standardized_log_residual", ascending=False
    )
    predictions.to_csv(OUTPUT_DIR / "regression_predictions.csv", index=False)
    predictions[predictions["inconsistency_flag"]].to_csv(
        OUTPUT_DIR / "regression_flagged_households.csv", index=False
    )
    pd.DataFrame(
        {"term": coefficient_names, "coefficient_log_welfare": full_coefficients}
    ).to_csv(OUTPUT_DIR / "regression_coefficients.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "survey": "MWI_1997",
                "records": len(data),
                "predictors": "urban + log1p(hsize) + subnatid",
                "folds": N_FOLDS,
                "weighted_r2_crossfit_log_welfare": weighted_r2(
                    target, crossfit_prediction_log, weights
                ),
                "weighted_rmse_crossfit_log_welfare": np.sqrt(
                    weighted_mean(residual_log**2, weights)
                ),
                "weighted_residual_sd_crossfit_log_welfare": residual_scale,
                "flag_threshold_standardized_abs_residual": FLAG_THRESHOLD,
                "flagged_records": int(data["inconsistency_flag"].sum()),
                "flagged_weight_share": weights[data["inconsistency_flag"]].sum()
                / weights.sum(),
            }
        ]
    )
    summary.to_csv(OUTPUT_DIR / "regression_summary.csv", index=False)

    plt.figure(figsize=(8, 6))
    flagged = data["inconsistency_flag"].to_numpy(bool)
    plt.scatter(
        target[~flagged],
        crossfit_prediction_log[~flagged],
        s=8,
        alpha=0.25,
        label="Not flagged",
    )
    plt.scatter(
        target[flagged],
        crossfit_prediction_log[flagged],
        s=14,
        color="#D55E00",
        alpha=0.75,
        label="Flagged",
    )
    limits = [min(target.min(), crossfit_prediction_log.min()), max(target.max(), crossfit_prediction_log.max())]
    plt.plot(limits, limits, color="black", linewidth=1, linestyle="--")
    plt.xlabel("Observed log welfare (2021 PPP USD/day)")
    plt.ylabel("Cross-fitted predicted log welfare")
    plt.title("Malawi 1997: observed versus predicted welfare")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "regression_observed_vs_predicted.png", dpi=180)
    plt.close()

    print(f"Valid records: {len(data)}")
    print(summary.to_string(index=False))
    print(f"Wrote regression outputs to {OUTPUT_DIR}")
    print(f"Coefficient count: {len(coefficient_names)}")


if __name__ == "__main__":
    main()