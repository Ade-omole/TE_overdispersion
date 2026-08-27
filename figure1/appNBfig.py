"""Plot mean and variance trajectories: approximate vs negative-binomial closure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----- Theory source
# "numeric"  -> integrated closure trajectories from the CSV (reproduces Fig. 1)
# "analytic" -> closed-form equilibria, drawn as flat reference lines
THEORY_SOURCE = "numeric"

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"
OUTPUT_PATH = FIG_DIR / "appNB.pdf"

# Input CSV produced by appNB.jl
CSV_BASENAME = "appNB_s:0.0001_N:10000_beta:0.026_delta:0.006_R:free_exp"

MAX_GENERATIONS_TO_PLOT = 2000  # None = all generations
SIM_BURN_IN = 2000
SIM_SAMPLE_EVERY = 100

csv_path = CSV_DIR / f"{CSV_BASENAME}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run appNB.jl first, or change CSV_BASENAME.")
    raise SystemExit(1)

data = pd.read_csv(csv_path)


# Rates used by the analytic expressions; must match the simulation CSV
U_RATE, V_RATE, S_SEL = 0.026, 0.006, 1e-4


def analytic_equilibria(u=U_RATE, v=V_RATE, s=S_SEL):
    """Closed-form equilibria from main_rev2 for Gaussian fitness w(n)=exp(-s n^2)."""
    mu_app = (u - v) * (1 - 2 * (u - v)) / (2 * s * (2 * u + 2 * v + 1))
    var_app = (u - v) / (2 * s)

    rho_nb = 1 + 3 * u + v
    mu_nb = (u - v) / (2 * s * rho_nb) - 0.5 - 3 * u - v
    var_nb = rho_nb * mu_nb

    skew_nb = (1 + 6 * u + 2 * v) / np.sqrt(var_nb)
    exkurt_nb = (1 + 6 * (3 * u + v) + 6 * (3 * u + v) ** 2) / var_nb

    return {
        "mu_app": mu_app, "var_app": var_app,
        "mu_nb": mu_nb, "var_nb": var_nb, "rho_nb": rho_nb,
        "skew_nb": skew_nb, "exkurt_nb": exkurt_nb,
        "charles": (u - v) / (2 * s),
    }


def theory_series(data, kind):
    """Return (mean, variance, charles) for one closure, honouring THEORY_SOURCE."""
    n = len(data)
    if THEORY_SOURCE == "analytic":
        eq = analytic_equilibria()
        m, var = (eq["mu_app"], eq["var_app"]) if kind == "approx" else (eq["mu_nb"], eq["var_nb"])
        const = lambda x: pd.Series(np.full(n, x), index=data.index)
        return const(m), const(var), const(eq["charles"])
    if kind == "approx":
        return data["Approx_Mean"], data["Approx_Variance"], data["Charles_Mean"]
    return data["NB_Mean"], data["NB_Variance"], data["Charles_Mean"]


mean_styles_approx = {
    "Simulated Mean": {"color": "#0173B2", "linestyle": "-", "linewidth": 5.5},
    "Approx Mean": {"color": "red", "linestyle": "-", "linewidth": 5.2},
}

variance_styles_approx = {
    "Simulated Variance": {"color": "#4B0082", "linestyle": "-", "linewidth": 1.9},
    "Approx Variance": {"color": "#000000", "linestyle": "-", "linewidth": 6.5},
    "Charles Mean": {"color": "#299438", "linestyle": "--", "linewidth": 6.5},
}

mean_styles_nb = {
    "Simulated Mean": {"color": "#0173B2", "linestyle": "-", "linewidth": 5.5},
    "NB Mean": {"color": "red", "linestyle": "-", "linewidth": 5.5},
}

variance_styles_nb = {
    "Simulated Variance": {"color": "#4B0082", "linestyle": "-", "linewidth": 1.8},
    "NB Variance": {"color": "#000000", "linestyle": "-", "linewidth": 6.5},
    "Charles Mean": {"color": "#299438", "linestyle": "--", "linewidth": 6.7},
}


def summarize_simulation(data):
    # Equilibrium values use the full run, not the plotting window:
    # samples every SIM_SAMPLE_EVERY generations after SIM_BURN_IN.
    idx = np.arange(SIM_BURN_IN, len(data), SIM_SAMPLE_EVERY)
    return {
        "Mean": data["Mean_X_Count"].iloc[idx].mean(),
        "Var.": data["Sim_Variance"].iloc[idx].mean(),
        "Skew.": data["Sim_Skewness"].iloc[idx].mean(),
        "Ex. Kurt.": data["Sim_Excess_kurtosis"].iloc[idx].mean(),
    }


def print_summary_table(data):
    sim_summary = summarize_simulation(data)
    if THEORY_SOURCE == "analytic":
        eq = analytic_equilibria()
        approx_mean, approx_variance = eq["mu_app"], eq["var_app"]
        nb_mean, nb_variance = eq["mu_nb"], eq["var_nb"]
        nb_skew, nb_ex_kurt = eq["skew_nb"], eq["exkurt_nb"]
        charles_mean = eq["charles"]
    else:
        approx_mean = data["Approx_Mean"].iloc[-1]
        approx_variance = data["Approx_Variance"].iloc[-1]
        nb_mean = data["NB_Mean"].iloc[-1]
        nb_variance = data["NB_Variance"].iloc[-1]
        nb_skew = data["NB_Skewness"].iloc[-1]
        nb_ex_kurt = data["NB_Excess_Kurtosis"].iloc[-1]
        charles_mean = data["Charles_Mean"].iloc[-1]

    summary_df = pd.DataFrame(
        {
            "Mean": [sim_summary["Mean"], nb_mean, approx_mean, charles_mean],
            "Var.": [sim_summary["Var."], nb_variance, approx_variance, charles_mean],
            "Skew.": [
                sim_summary["Skew."],
                nb_skew,
                np.nan,
                1 / np.sqrt(charles_mean) if charles_mean > 0 else np.nan,
            ],
            "Ex. Kurt.": [
                sim_summary["Ex. Kurt."],
                nb_ex_kurt,
                np.nan,
                1 / charles_mean if charles_mean > 0 else np.nan,
            ],
        },
        index=[
            "Stochastic sim.",
            "Theory (neg. bin.)",
            "Theory (approx.)",
            "Ch. & Ch. (1983)",
        ],
    )

    print("\nSummary table:\n")
    print(summary_df.to_string(float_format=lambda x: f"{x:.2f}", na_rep="---"))


def add_legend(ax):
    ax.legend(loc="lower right", prop={"size": 13, "weight": "bold"}, framealpha=0.9)


def panel_label(ax, letter):
    ax.text(
        -0.06, 1.06, letter, transform=ax.transAxes,
        fontsize=24, fontweight="bold", va="top", ha="right",
    )


def plot_approx(data, ax, title):
    if MAX_GENERATIONS_TO_PLOT is not None:
        data = data.iloc[:MAX_GENERATIONS_TO_PLOT]

    generation = data["generation"]
    simulated_mean = data["Mean_X_Count"]
    simulated_variance = data["Sim_Variance"]
    approx_mean, approx_variance, charles_mean = theory_series(data, "approx")

    ax.plot(generation, approx_mean, label="Theoretical Mean", **mean_styles_approx["Approx Mean"])
    ax.plot(
        generation,
        simulated_mean,
        label="Simulated Mean",
        **mean_styles_approx["Simulated Mean"],
        alpha=0.9,
    )
    ax.plot(
        generation,
        simulated_variance,
        label="Simulated Variance",
        **variance_styles_approx["Simulated Variance"],
        alpha=0.85,
    )
    ax.plot(
        generation,
        approx_variance,
        label="Theoretical Variance",
        **variance_styles_approx["Approx Variance"],
    )
    ax.plot(
        generation,
        charles_mean,
        label="Charlesworth (1983) Mean",
        **variance_styles_approx["Charles Mean"],
    )

    ax.set_xlabel("Generation", fontsize=22, fontweight="bold")
    ax.set_ylabel("Copy number", fontsize=25, fontweight="bold")
    ax.set_title(title, fontsize=20, fontweight="bold", pad=5)
    ax.grid(True, linestyle="--", alpha=0, linewidth=0.5)
    ax.tick_params(axis="both", which="major", labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    add_legend(ax)
    panel_label(ax, "a")


def plot_nb(data, ax, title):
    if MAX_GENERATIONS_TO_PLOT is not None:
        data = data.iloc[:MAX_GENERATIONS_TO_PLOT]

    generation = data["generation"]
    simulated_mean = data["Mean_X_Count"]
    simulated_variance = data["Sim_Variance"]
    nb_mean, nb_variance, charles_mean = theory_series(data, "nb")

    ax.plot(generation, nb_mean, label="Theoretical Mean", **mean_styles_nb["NB Mean"])
    ax.plot(
        generation,
        simulated_mean,
        label="Simulated Mean",
        **mean_styles_nb["Simulated Mean"],
        alpha=0.9,
    )
    ax.plot(
        generation,
        simulated_variance,
        label="Simulated Variance",
        **variance_styles_nb["Simulated Variance"],
        alpha=0.85,
    )
    ax.plot(
        generation,
        nb_variance,
        label="Theoretical Variance",
        **variance_styles_nb["NB Variance"],
    )
    ax.plot(
        generation,
        charles_mean,
        label="Charlesworth (1983) Mean",
        **variance_styles_nb["Charles Mean"],
    )

    ax.set_xlabel("Generation", fontsize=22, fontweight="bold")
    ax.set_title(title, fontsize=20, fontweight="bold", pad=5)
    ax.grid(True, linestyle="--", alpha=0, linewidth=0.5)
    ax.tick_params(axis="both", which="major", labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    add_legend(ax)
    panel_label(ax, "b")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7), dpi=120)

plot_approx(data, ax1, r"$\bf{Approximate\ closure}$")
plot_nb(data, ax2, r"$\bf{Negative\ binomial\ closure}$")
print_summary_table(data)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
print(f"\nSaved figure to: {OUTPUT_PATH}")
plt.show()
