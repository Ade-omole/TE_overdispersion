"""Plot dispersion ratio rho versus map length R (varying-R sweep)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, LogFormatterMathtext, NullLocator
from scipy.interpolate import make_interp_spline

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"
OUTPUT_PATH = FIG_DIR / "varyingRFig.pdf"

# Input CSV basenames (without fitness suffix); produced by varyingRsim.jl.
CSV_BASENAME_PRIMARY = "sweep_R_roze_topright_burnin_binomial"
CSV_BASENAME_SIM_ECT = "sweep_R_ectFitRoze_burnin_binomial"
FITNESS = "exp"  # "exp" or "syn"

csv_path = CSV_DIR / f"{CSV_BASENAME_PRIMARY}_{FITNESS}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run varyingRsim.jl first, or change CSV_BASENAME_PRIMARY / FITNESS.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

csv_path_ect = CSV_DIR / f"{CSV_BASENAME_SIM_ECT}_{FITNESS}.csv"
if not csv_path_ect.is_file():
    print(f"CSV not found: {csv_path_ect}")
    print("Add the ectFitRoze sweep CSV or change CSV_BASENAME_SIM_ECT / FITNESS.")
    raise SystemExit(1)

df_ect = pd.read_csv(csv_path_ect)
print(f"Loaded {csv_path_ect}, shape {df_ect.shape}")
if len(df["R"]) != len(df_ect["R"]) or not np.allclose(df["R"].values, df_ect["R"].values):
    print("Warning: R grids differ between primary and ectFitRoze CSVs; ect simulation uses df_ect R/sim_p.")

R = df["R"].values

sim_p = df["sim_p"].values
sim_p_ect = df_ect["sim_p"].values
app_p = df["app_p"].values
nb_p = df["nb_p"].values
has_rho_numerical = "rho_numerical" in df.columns
if has_rho_numerical:
    rho_numerical = df["rho_numerical"].values

# Plot styling
LEFT_SIM_MARKER_SIZE = 7
THEORY_LINEWIDTH = 4.0
MARKER_EDGEWIDTH = 2.5

TICK_LABELSIZE = 15
LEGEND_FONTSIZE = 16
AXIS_LABELSIZE = 20

SHOW_VERTICAL_GRID = True
VERTICAL_GRID_ALPHA = 0.18
VERTICAL_GRID_LINESTYLE = "-"
VERTICAL_GRID_LINEWIDTH = 0.8

SHOW_HORIZONTAL_GRID = True
HORIZONTAL_GRID_ALPHA = 0.18
HORIZONTAL_GRID_LINESTYLE = "-"
HORIZONTAL_GRID_LINEWIDTH = 0.8
Y_MINOR_SUBDIVISIONS = 5
HORIZONTAL_MINOR_GRID_ALPHA = 0.10
HORIZONTAL_MINOR_GRID_LINESTYLE = ":"
HORIZONTAL_MINOR_GRID_LINEWIDTH = 0.55

idx_sort = np.argsort(R)
x_line = R[idx_sort]


def smooth_log_spline(x_sorted, y_sorted, n_interp=500):
    """
    Fit a cubic spline in log10(x) space and evaluate on a fine grid,
    returning (x_fine, y_fine) for smooth plotting on a log x-axis.
    Requires at least 4 unique points; falls back to raw data if fewer.
    """
    log_x = np.log10(x_sorted)
    _, unique_idx = np.unique(log_x, return_index=True)
    log_x_u = log_x[unique_idx]
    y_u = y_sorted[unique_idx]
    if len(log_x_u) < 4:
        return x_sorted, y_sorted
    spline = make_interp_spline(log_x_u, y_u, k=3)
    log_x_fine = np.linspace(log_x_u[0], log_x_u[-1], n_interp)
    y_fine = spline(log_x_fine)
    return 10**log_x_fine, y_fine


_log_x_major_fallback = LogFormatterMathtext(base=10)


def _format_log_x_major(x, pos):
    """Log x-axis: show 1 instead of 10^0; other ticks use matplotlib's default log math text."""
    if x <= 0 or not np.isfinite(x):
        return ""
    if np.isclose(x, 1.0, rtol=1e-12, atol=0.0):
        return "1"
    return _log_x_major_fallback(x, pos)


color_nb = "#2ca02c"
color_app = "#9467bd"
color_sim_p = "#d62728"
color_sim_ect = "#B8860B"
color_roze = "#0173B2"

line_kwargs = dict(
    linestyle="-",
    linewidth=THEORY_LINEWIDTH,
    zorder=3,
    antialiased=True,
    solid_joinstyle="round",
    solid_capstyle="round",
)

fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

inv_app_p = np.array([1.0 / max(p, 1e-12) for p in app_p])
inv_nb_p = np.array([1.0 / max(p, 1e-12) for p in nb_p])

x_smooth_app, y_smooth_app = smooth_log_spline(x_line, inv_app_p[idx_sort])
x_smooth_nb, y_smooth_nb = smooth_log_spline(x_line, inv_nb_p[idx_sort])
ax.plot(x_smooth_app, y_smooth_app, color=color_app, label="Approximate", **line_kwargs)
ax.plot(x_smooth_nb, y_smooth_nb, color=color_nb, label="Parametric", **line_kwargs)

if has_rho_numerical:
    x_roze_s, y_roze_s = smooth_log_spline(x_line, rho_numerical[idx_sort])
    ax.plot(x_roze_s, y_roze_s, color=color_roze, label="Roze (2023)", **line_kwargs)

ax.scatter(
    R,
    1.0 / sim_p,
    marker="o",
    s=LEFT_SIM_MARKER_SIZE**2,
    zorder=12,
    color=color_sim_p,
    edgecolors=color_sim_p,
    linewidths=MARKER_EDGEWIDTH,
    label=r"Simulation ($w = e^{-sn^2}$)",
)
R_ect = df_ect["R"].values
ax.scatter(
    R_ect,
    1.0 / sim_p_ect,
    marker="o",
    s=LEFT_SIM_MARKER_SIZE**2,
    zorder=10,
    color=color_sim_ect,
    edgecolors=color_sim_ect,
    linewidths=MARKER_EDGEWIDTH,
    label=r"Simulation ($w = e^{-\beta n_p}$)",
)

ax.set_xscale("log")
ax.xaxis.set_major_formatter(FuncFormatter(_format_log_x_major))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlabel(r"map length (R)", fontsize=AXIS_LABELSIZE, fontweight="bold")
ax.tick_params(axis="x", labelsize=TICK_LABELSIZE)
fig.canvas.draw()
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
ax.set_ylabel(r"Dispersion ratio ($\rho$)", fontsize=AXIS_LABELSIZE, fontweight="bold")
ax.tick_params(axis="y", labelsize=TICK_LABELSIZE)
for label in ax.get_yticklabels():
    label.set_fontweight("bold")

ax.yaxis.set_minor_locator(AutoMinorLocator(Y_MINOR_SUBDIVISIONS))
ax.tick_params(axis="y", which="minor", length=4, width=0.7)

ax.grid(
    SHOW_HORIZONTAL_GRID,
    which="major",
    axis="y",
    linestyle=HORIZONTAL_GRID_LINESTYLE,
    linewidth=HORIZONTAL_GRID_LINEWIDTH,
    alpha=HORIZONTAL_GRID_ALPHA,
)
ax.grid(
    SHOW_HORIZONTAL_GRID,
    which="minor",
    axis="y",
    linestyle=HORIZONTAL_MINOR_GRID_LINESTYLE,
    linewidth=HORIZONTAL_MINOR_GRID_LINEWIDTH,
    alpha=HORIZONTAL_MINOR_GRID_ALPHA,
)
ax.grid(
    SHOW_VERTICAL_GRID,
    which="major",
    axis="x",
    linestyle=VERTICAL_GRID_LINESTYLE,
    linewidth=VERTICAL_GRID_LINEWIDTH,
    alpha=VERTICAL_GRID_ALPHA,
)
ax.legend(loc="best", prop={"size": LEGEND_FONTSIZE, "weight": "bold"})

plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
print(f"Saved figure: {OUTPUT_PATH}")
plt.show()
