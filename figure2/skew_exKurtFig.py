"""Plot simulated vs NB-theory skewness and excess kurtosis over generations."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"
OUTPUT_PATH = FIG_DIR / "skew_and_exkurt.pdf"

# Input CSV produced by appNB.jl (same as appNBfig.py)
CSV_BASENAME = "appNB_s:0.0001_N:10000_beta:0.026_delta:0.006_R:free_exp"

MAX_GENERATIONS_TO_PLOT = 1000  # None = all generations
SIM_BURN_IN = 2000
SIM_SAMPLE_EVERY = 100

csv_path = CSV_DIR / f"{CSV_BASENAME}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run appNB.jl first, or change CSV_BASENAME.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

plt.style.use("default")
plt.rcParams.update(
    {
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "legend.fontsize": 12,
    }
)

if MAX_GENERATIONS_TO_PLOT is not None:
    df_plot = df.iloc[:MAX_GENERATIONS_TO_PLOT].reset_index(drop=True)
else:
    df_plot = df.reset_index(drop=True)

x = np.arange(len(df_plot))

colors = ["#1f77b4", "#ff7f0e", "#000000"]
linestyles = ["-", "-", "-"]
linewidths = [2.5, 6.5, 5.5]

label_map = {
    "Sim_Skewness": "Simulated Skewness",
    "NB_Skewness": "Theoretical (NB) Skewness",
    "Sim_Excess_kurtosis": "Simulated Excess Kurtosis",
    "NB_Excess_Kurtosis": "Theoretical (NB) Excess Kurtosis",
}


def plot_data_group(data, ax, columns, title, ylabel):
    for i, col in enumerate(columns):
        if col in data.columns:
            ax.plot(
                x,
                data[col],
                color=colors[i % len(colors)],
                linestyle=linestyles[i % len(linestyles)],
                linewidth=linewidths[i % len(linewidths)],
                alpha=0.9,
                label=label_map.get(col, col),
                marker="o" if len(data) < 50 else None,
                markersize=4 if len(data) < 50 else 0,
            )
        else:
            print(f"Warning: column '{col}' not found in the dataset")

    ax.set_xlabel("Generation", fontsize=22, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=22, fontweight="bold")
    ax.set_title(title, fontsize=20, fontweight="bold", pad=5)
    ax.grid(True, alpha=0)
    ax.legend(loc="best", prop={"size": 14, "weight": "bold"})
    ax.tick_params(axis="both", which="major", labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

    print(f"\n{title}:")
    for col in columns:
        if col in data.columns:
            print(f"  {label_map.get(col, col)}: {data[col].iloc[-1]}")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

plot_data_group(
    df_plot,
    ax1,
    ["Sim_Skewness", "NB_Skewness"],
    "Skewness (NB)",
    "Skewness Value",
)
plot_data_group(
    df_plot,
    ax2,
    ["Sim_Excess_kurtosis", "NB_Excess_Kurtosis"],
    "Excess Kurtosis (NB)",
    "Excess Kurtosis Value",
)

# Piecewise averages for simulated skewness and excess kurtosis (same burn-in as appNBfig.py)
if MAX_GENERATIONS_TO_PLOT is not None:
    df_summary = df.iloc[:MAX_GENERATIONS_TO_PLOT]
else:
    df_summary = df

skewness = df_summary["Sim_Skewness"]
excess_kurtosis = df_summary["Sim_Excess_kurtosis"]

indices_skew = np.arange(SIM_BURN_IN, len(skewness), SIM_SAMPLE_EVERY)
indices_kurt = np.arange(SIM_BURN_IN, len(excess_kurtosis), SIM_SAMPLE_EVERY)

final_average_sim_skewness = skewness.iloc[indices_skew].mean()
final_average_sim_kurtosis = excess_kurtosis.iloc[indices_kurt].mean()

print(f"\nSimulated skewness (piecewise average): {final_average_sim_skewness:.2f}")
print(f"Simulated excess kurtosis (piecewise average): {final_average_sim_kurtosis:.2f}")

plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
print(f"Figure saved to: {OUTPUT_PATH}")
plt.show()
