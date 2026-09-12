"""Fit weighted GB2 models to the Malawi 1997 upper tail.

The script fits the GB2 distribution to the full positive welfare sample and
to samples after removing the largest k observations. It then compares the
weighted empirical survival function with the fitted GB2 survival function at
the largest observed values.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, betaincc, digamma


ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = ROOT / "01-input" / "country" / "MWI_1997.dta"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"

TOP_K_VALUES = (0, 1, 2, 5, 10, 20)
TAIL_COUNT = 30
STARTS = ((1.0, 1.0, 2.0), (1.0, 2.0, 3.0), (2.0, 1.0, 2.0), (0.5, 1.0, 1.5), (2.0, 2.0, 5.0))


def load_malawi_data():
    data = pd.read_stata(
        INPUT_PATH,
        columns=["welfare", "weight", "cpi2021", "icp2021"],
    )
    welfare = data["welfare"] / (data["cpi2021"] * data["icp2021"] * 365)
    weight = data["weight"]
    valid = (
        np.isfinite(welfare)
        & welfare.gt(0)
        & np.isfinite(weight)
        & weight.gt(0)
    )
    return welfare.loc[valid].to_numpy(float), weight.loc[valid].to_numpy(float)


def weighted_quantile(values, weights, probability):
    order = np.argsort(values)
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order]) / np.sum(weights)
    return sorted_values[np.searchsorted(cumulative, probability, side="left")]


def gb2_logpdf(log_values, log_a, log_b, log_p, log_q):
    a, b, p, q = np.exp([log_a, log_b, log_p, log_q])
    centered = log_values - log_b
    return (
        log_a
        - log_b
        - betaln(p, q)
        + (a * p - 1) * centered
        - (p + q) * np.logaddexp(0.0, a * centered)
    )


def fit_gb2(values, weights):
    log_values = np.log(values)
    weighted_median = weighted_quantile(values, weights, 0.5)
    normalized_weights = weights / np.sum(weights)

    def objective(log_parameters):
        log_density = gb2_logpdf(log_values, *log_parameters)
        value = -np.sum(normalized_weights * log_density)
        return value if np.isfinite(value) else 1e100

    def objective_with_gradient(log_parameters):
        log_a, log_b, log_p, log_q = log_parameters
        a, b, p, q = np.exp(log_parameters)
        centered = log_values - log_b
        scaled = a * centered
        logistic = 1.0 / (1.0 + np.exp(-np.clip(scaled, -700, 700)))
        log_density = gb2_logpdf(log_values, log_a, log_b, log_p, log_q)
        value = -np.sum(normalized_weights * log_density)
        if not np.isfinite(value):
            return 1e100, np.zeros(4)

        # Derivatives with respect to log(a), log(b), log(p), and log(q).
        log_one_plus = np.logaddexp(0.0, scaled)
        gradient = np.array(
            [
                1 + a * p * centered - (p + q) * a * centered * logistic,
                -a * p + (p + q) * a * logistic,
                p * (a * centered - digamma(p) + digamma(p + q) - log_one_plus),
                q * (-digamma(q) + digamma(p + q) - log_one_plus),
            ]
        )
        # The expression above is per-observation; apply the negative
        # log-likelihood sign and sampling weights.
        weighted_gradient = -np.sum(normalized_weights[None, :] * gradient, axis=1)
        return value, weighted_gradient

    best = None
    for a_start, p_start, q_start in STARTS:
        initial = np.log([a_start, weighted_median, p_start, q_start])
        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=lambda parameters: objective_with_gradient(parameters)[1],
            bounds=[(-6, 6), (np.log(1e-4), np.log(1e6)), (-6, 6), (-6, 6)],
            options={"maxiter": 500, "ftol": 1e-10},
        )
        if result.success and (best is None or result.fun < best.fun):
            best = result

    if best is None:
        raise RuntimeError("GB2 optimization failed for all starting values")

    a, b, p, q = np.exp(best.x)
    return {
        "a": a,
        "b": b,
        "p": p,
        "q": q,
        "weighted_neg_loglik": best.fun,
        "converged": bool(best.success),
        "optimizer_message": str(best.message),
    }


def gb2_survival(values, parameters):
    a, b, p, q = (parameters[name] for name in ("a", "b", "p", "q"))
    transformed = (values / b) ** a
    z = transformed / (1 + transformed)
    return betaincc(p, q, z)


def weighted_mean(values, weights):
    return np.sum(values * weights) / np.sum(weights)


def main():
    values, weights = load_malawi_data()
    order = np.argsort(values)[::-1]
    rows = []
    tail_rows = []
    fits = {}

    for top_k in TOP_K_VALUES:
        keep = np.ones(len(values), dtype=bool)
        if top_k:
            keep[order[:top_k]] = False
        fit_values = values[keep]
        fit_weights = weights[keep]
        parameters = fit_gb2(fit_values, fit_weights)
        fits[top_k] = parameters
        rows.append(
            {
                "survey": "MWI_1997",
                "removed_top_records": top_k,
                "fit_records": len(fit_values),
                "fit_weighted_mean_ppp_day": weighted_mean(fit_values, fit_weights),
                **parameters,
            }
        )

        tail_order = order[:TAIL_COUNT]
        empirical_survival = np.array(
            [np.sum(weights[values >= value]) / np.sum(weights) for value in values[tail_order]]
        )
        fitted_survival = gb2_survival(values[tail_order], parameters)
        for rank, index in enumerate(tail_order, start=1):
            tail_rows.append(
                {
                    "survey": "MWI_1997",
                    "removed_top_records": top_k,
                    "rank_largest_first": rank,
                    "observed_welfare_ppp_day": values[index],
                    "sampling_weight": weights[index],
                    "empirical_survival_full_sample": empirical_survival[rank - 1],
                    "gb2_survival_at_observed_value": fitted_survival[rank - 1],
                    "empirical_to_gb2_survival_ratio": empirical_survival[rank - 1] / fitted_survival[rank - 1]
                    if fitted_survival[rank - 1] > 0
                    else np.nan,
                }
            )

    summary = pd.DataFrame(rows)
    tail = pd.DataFrame(tail_rows)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_DIR / "gb2_fit_sensitivity.csv", index=False)
    tail.to_csv(OUTPUT_DIR / "gb2_observed_vs_fitted_tail.csv", index=False)

    full_parameters = fits[0]
    grid = np.geomspace(values.min(), values.max(), 400)
    empirical = np.array([np.sum(weights[values >= point]) / np.sum(weights) for point in grid])
    plt.figure(figsize=(9, 6))
    plt.step(grid, empirical, where="post", color="#0072B2", label="Observed weighted survival")
    for top_k, parameters in fits.items():
        label = "GB2: full sample" if top_k == 0 else f"GB2: excluding top {top_k}"
        plt.plot(grid, gb2_survival(grid, parameters), label=label)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Welfare (2021 PPP USD/day, log scale)")
    plt.ylabel("Survival probability (log scale)")
    plt.title("Malawi 1997: observed upper tail versus GB2 fits")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "gb2_observed_vs_fitted_tail.png", dpi=180)
    plt.close()

    print(f"Valid records: {len(values)}")
    print(summary[["removed_top_records", "a", "b", "p", "q", "weighted_neg_loglik"]].to_string(index=False))
    print(f"Wrote GB2 outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()