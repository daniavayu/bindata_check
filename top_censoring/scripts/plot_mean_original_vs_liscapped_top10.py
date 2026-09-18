"""Create a dumbbell plot for the ten surveys with the largest mean change."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT_PATH = ROOT / "top_censoring" / "outputs" / "all_surveys_four_indicators_lis_comparison.csv"
OUTPUT_PATH = ROOT / "top_censoring" / "outputs" / "plot_mean_original_vs_liscapped_top10_change.png"


def main():
    data = pd.read_csv(INPUT_PATH)
    change_column = "mean_ppp_day change (%)"
    data = data.loc[data[change_column].notna()].copy()
    data["label"] = (
        data["Country"].astype(str)
        + "-"
        + data["Year"].astype(int).astype(str)
        + " ("
        + data["Survey"].astype(str)
        + ")"
    )

    # The capping changes are negative; the most influential surveys are the
    # observations with the largest absolute reduction in the mean.
    top = data.nsmallest(10, change_column).sort_values(change_column)

    figure, axis = plt.subplots(figsize=(12, 7.5))
    positions = range(len(top))
    original = top["mean_ppp_day original"]
    capped = top["mean_ppp_day LIS-capped"]

    for position, original_value, capped_value in zip(positions, original, capped):
        axis.plot([capped_value, original_value], [position, position], color="#777777", linewidth=1.2, zorder=1)

    axis.scatter(original, positions, color="#1f77b4", s=28, label="Original mean", zorder=2)
    axis.scatter(capped, positions, color="#d62728", s=28, label="LIS-capped mean", zorder=3)

    for position, value in zip(positions, top[change_column]):
        x = max(original.iloc[position], capped.iloc[position])
        axis.text(
            x + 0.55,
            position,
            f"{abs(value):.1f}%",
            va="center",
            fontsize=10,
            color="#444444",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8},
        )

    axis.set_yticks(list(positions))
    axis.set_yticklabels(top["label"])
    axis.set_xlabel("Mean PPP/day")
    axis.set_title("Top 10 surveys by largest % change in mean (Original vs LIS-capped)")
    axis.grid(axis="x", color="#dddddd", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.legend(loc="lower right", frameon=True)
    figure.tight_layout()
    figure.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    print(f"Wrote {OUTPUT_PATH}")
    print(top[["Country", "Year", "Survey", "mean_ppp_day original", "mean_ppp_day LIS-capped", change_column]].to_string(index=False))


if __name__ == "__main__":
    main()