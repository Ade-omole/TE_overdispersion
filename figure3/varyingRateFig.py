"""Plot dispersion ratio and higher moments versus transposition rate u."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from scipy import integrate
from scipy.interpolate import make_interp_spline

# Theory source for the NB curves
NB_THEORY_SOURCE = "analytic"  # "numeric" or "analytic"

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"
OUTPUT_PATH = FIG_DIR / "varyingRateFig.pdf"

R_val = 10.0
# Roze pre-transposition curve
ROZE_PRETRANSP_MODE = "simple"  # "simple" or "recomb"
fitness = "exp"
EXCLUDE_LAST_N = 3  # e.g. set to 3 to drop the last 3 data points
csv_path = CSV_DIR / "sweep_u_more_values_vConst:0.0001_6_binomial_exp.csv" ## R2
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Make sure the CSV exists at that exact path.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
if EXCLUDE_LAST_N > 0:
    df = df.iloc[:-EXCLUDE_LAST_N].copy()
print(f"Loaded {csv_path}, shape {df.shape}")
v_label = df["v"].iloc[0] if "v" in df.columns else "0.00001"

# Post-transposition columns
sim_mean = df["sim_μ"].values
app_mean = df["app_μ"].values
nb_mean = df["nb_μ"].values
sim_var = df["sim_σ²"].values
app_var = df["app_σ²"].values
nb_var = df["nb_σ²"].values
sim_p = df["sim_p"].values
app_p = df["app_p"].values
nb_p = df["nb_p"].values
has_skew = "sim_skew" in df.columns
if has_skew:
    sim_skew = df["sim_skew"].values
    nb_skew = df["nb_skew"].values
    sim_exkurt = df["sim_exkurt"].values
    nb_exkurt = df["nb_exkurt"].values

# Pre-transposition (m0) columns
has_m0_sim = "sim_p_m0" in df.columns
if has_m0_sim:
    sim_p_m0 = df["sim_p_m0"].values

# Style
LEFT_SIM_MARKER_SIZE = 8  # Left-panel simulation markers
RIGHT_SIM_MARKER_SIZE = 11  # Right-panel simulation markers

MAIN_LINEWIDTH = 5  # Left-panel theory curves: approximate, parametric, and Roze
SECONDARY_LINEWIDTH = 5.0  # Right-panel parametric curves: skewness and excess kurtosis
MARKER_EDGEWIDTH = 3.5  # Outline thickness of the simulation point markers

X_TICK_LABELSIZE = 15  # X-axis tick-label size for both panels
Y_TICK_LABELSIZE = 15  # Y-axis tick-label size for both panels

LEGEND_FONTSIZE = 16  # Legend text size for the right panel
LEFT_LEGEND_FONTSIZE = 14  # Left-panel legend text size; slightly smaller since labels are longer

X_AXIS_LABELSIZE = 18  # X-axis label size for both panels
Y_AXIS_LABELSIZE = 20  # Y-axis label size for both panels

# Vertical grid
SHOW_VERTICAL_GRID = True
VERTICAL_GRID_ALPHA = 0.18
VERTICAL_GRID_LINESTYLE = "-"
VERTICAL_GRID_LINEWIDTH = 0.8

# Horizontal grid
SHOW_HORIZONTAL_GRID = True
HORIZONTAL_GRID_ALPHA = 0.18
HORIZONTAL_GRID_LINESTYLE = "-"
HORIZONTAL_GRID_LINEWIDTH = 0.8


def compute_E1(R, u, eps=1e-12):
    """E₁ (ε₁) for linear genetic map of length R Morgans (Roze)."""
    R = max(R, eps)
    rho_R = u / R
    if R <= 1.0:
        half_plus_rho = 0.5 + rho_R
        if rho_R <= 0.0 or half_plus_rho <= 0.0:
            return 0.0
        return (2.0 / R) * (
            (1.0 + 2.0 * rho_R) * (np.log(half_plus_rho) - np.log(max(rho_R, eps))) - 1.0
        )

    def integrand(x):
        denom = (1.0 - np.exp(-2.0 * x)) / 2.0 + 2.0 * u
        return (R - x) / max(np.finfo(float).eps, denom)

    res, _ = integrate.quad(integrand, 0.0, R, limit=200)
    return (2.0 / (R**2)) * res


def compute_rho(E1, u, v, alpha=0.0, eps=1e-10):
    """ρ from Roze: ρ ≈ 1 + (ε₁/(1-uε₁)) * ((u+v+α)/2)."""
    denom = 1.0 - u * E1
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    return 1.0 + (E1 / denom) * ((u + v + alpha) / 2.0)


S_VAL = 1e-4  # Gaussian/exponential fitness w(n)=exp(-s n^2), matches fit="exp" and the simulation


def nb_analytic(u, v, s):
    """Closed-form negative-binomial equilibrium (main_rev2, Gaussian fitness)."""
    rho = 1.0 + 3.0 * u + v + s                                  # Eq. (14)
    rho15 = 1.0 + 3.0 * u + v                                    # rho used in Eq. (15)
    mu = (u - v) / (2.0 * s * rho15) - 0.5 - 3.0 * u - v         # Eq. (15)
    var = rho15 * mu                                             # Eq. (15)
    skew = (1.0 + 6.0 * u + 2.0 * v) / np.sqrt(var)              # Eq. (16)
    exkurt = (1.0 + 6.0 * (3.0 * u + v) + 6.0 * (3.0 * u + v) ** 2) / var  # Eq. (16)
    # Eq. (23), unexpanded ratio
    rho_pre = (mu + var * (1.0 - 3.0 * (u - v)) + s * var**2 / mu) / (
        2.0 * mu * (1.0 - (u - v))
    )
    return rho, mu, var, skew, exkurt, rho_pre


def rho_m0(mu, sig2, s=S_VAL):
    """Pre-transposition dispersion ratio rho_{m0} = Var[m0]/E[m0] (Eqs. 6-7)."""
    p = mu / sig2
    sigma = np.sqrt(sig2)
    gamma = ((2.0 - p) / p) / sigma          # Eq. 9-10 (theta/sigma)
    kappa = (6.0 * (1.0 - p) + p**2) / sig2 + 3.0  # Eq. 11 (alpha/sigma^2 + 3)

    beta1, beta2 = -2.0 * s * mu, -2.0 * s
    E_m0 = mu + beta1 * sig2 + 0.5 * beta2 * gamma * sig2**1.5
    Var_m0 = (
        0.5 * (mu + sig2)
        + 0.5 * beta1 * (sig2 + gamma * sig2**1.5)
        + 0.25 * beta2 * (gamma * sig2**1.5 + (kappa - 1.0) * sig2**2)
    )
    return Var_m0 / E_m0


u_vals = df["u"].values
v_vals = df["v"].values
E1_vals = np.array([compute_E1(R_val, u) for u in u_vals])

# Recombination-map based Roze prediction (unused by default)
rho_roze_recomb = np.array([compute_rho(E1, u, v) for E1, u, v in zip(E1_vals, u_vals, v_vals)])

# Roze pre-transposition: rho = 1+u
rho_roze_simple = 1.0 + u_vals + v_vals

rho_roze = rho_roze_recomb if ROZE_PRETRANSP_MODE == "recomb" else rho_roze_simple

# Roze (2023), post-transposition: rho = 1 + 3u
rho_roze_3u = 1.0 + 3.0 * u_vals 

# Our pre-transposition prediction
rho_m0_vals = np.array([rho_m0(mu, sig2) for mu, sig2 in zip(nb_mean, nb_var)])

if NB_THEORY_SOURCE == "analytic":
    # Override every negative-binomial quantity with its closed form; the
    _an = np.array([nb_analytic(u, v, S_VAL) for u, v in zip(u_vals, v_vals)])
    _rho, nb_mean, nb_var, _sk, _ek, rho_m0_vals = (_an[:, k] for k in range(6))
    nb_p = 1.0 / _rho
    if has_skew:
        nb_skew, nb_exkurt = _sk, _ek

idx_sort_u = np.argsort(u_vals)
x_line_u = u_vals[idx_sort_u]

# ----- Print Julia-style table per (u, v)
n_total = len(u_vals)
for i in range(n_total):
    u, v = u_vals[i], v_vals[i]
    print(f"\n[{i+1}/{n_total}]  u = {u}  (v = {v}, fitness = {fitness}) [Roze]")
    print("  ┌─────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐")
    print("  │             │     Mean     │   Variance   │      p       │   Skewness   │  Ex. Kurtosis│")
    print("  ├─────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤")
    sm, sv, sp = sim_mean[i], sim_var[i], sim_p[i]
    if has_skew:
        print(f"  │ Simulation  │ {sm:12.2f} │ {sv:12.2f} │ {sp:12.4f} │ {sim_skew[i]:12.2f} │ {sim_exkurt[i]:12.2f} │")
    else:
        print(f"  │ Simulation  │ {sm:12.2f} │ {sv:12.2f} │ {sp:12.4f} │      —       │      —       │")

    am, av, ap = app_mean[i], app_var[i], app_p[i]
    print(f"  │ Approx      │ {am:12.2f} │ {av:12.2f} │ {ap:12.4f} │      —       │      —       │")

    nm, nv, np_ = nb_mean[i], nb_var[i], nb_p[i]
    if has_skew:
        print(f"  │ Theory (NB) │ {nm:12.2f} │ {nv:12.2f} │ {np_:12.4f} │ {nb_skew[i]:12.2f} │ {nb_exkurt[i]:12.2f} │")
    else:
        print(f"  │ Theory (NB) │ {nm:12.2f} │ {nv:12.2f} │ {np_:12.4f} │      —       │      —       │")

    e1, rho = E1_vals[i], rho_roze[i]
    print(f"  │ Roze        │ E₁={e1:.4f}   │ ρ={rho:.4f}     │      —       │      —       │      —       │")
    print("  └─────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘")
    print(
        f"  ρ: sim = {1.0/max(sim_p[i], 1e-12):.4f},  "
        f"approx = {1.0/max(app_p[i], 1e-12):.4f},  "
        f"parametric = {1.0/max(nb_p[i], 1e-12):.4f},  "
        f"Roze = {rho_roze[i]:.4f}"
    )


def smooth_log_spline(x_sorted, y_sorted, n_interp=500):
    """Fit a cubic spline in log10(x) space for smooth plotting."""
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


def _format_log_tick(value, _pos):
    exponent = int(np.floor(np.log10(value)))
    coeff = value / (10**exponent)
    coeff_rounded = int(round(coeff))
    if np.isclose(coeff, 1.0):
        return rf"$\mathbf{{10^{{{exponent}}}}}$"
    return rf"$\mathbf{{{coeff_rounded} \times 10^{{{exponent}}}}}$"


def style_log_xaxis(ax, x_values):
    ax.set_xscale("log")

    candidate_ticks = np.array([1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    x_min = np.min(x_values)
    x_max = np.max(x_values)
    visible_ticks = candidate_ticks[(candidate_ticks >= x_min) & (candidate_ticks <= x_max)]
    if len(visible_ticks) == 0:
        visible_ticks = candidate_ticks

    ax.xaxis.set_major_locator(FixedLocator(visible_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(_format_log_tick))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", which="major", labelsize=X_TICK_LABELSIZE, length=5, width=1.0, direction="out")


def plot_left_panel(ax):
    ax.plot(x_smooth_app, y_smooth_app, color=color_app, linestyle="--", label="Approximate (Post-transp.)", **line_kwargs)
    ax.plot(x_smooth_nb, y_smooth_nb, color=color_nb, linestyle="-", label="Parametric (Post-transp.)", **line_kwargs)
    ax.plot(
        x_smooth_m0,
        y_smooth_m0,
        color=color_m0,
        linestyle="-",
        label=r"Parametric (Pre-transp.)",
        **line_kwargs,
    )
    # Drawn after (and with a higher zorder than) the Parametric (Pre-transp.) line
    ax.plot(
        x_smooth_roze,
        y_smooth_roze,
        color=color_roze,
        linestyle=":",
        label="Roze (2023)-Pre-transp.",
        **{**line_kwargs, "zorder": line_kwargs["zorder"] + 1},
    )
    ax.plot(
        x_smooth_roze3u,
        y_smooth_roze3u,
        color=color_roze3u,
        linestyle=":",
        label="Roze (2023)-Post-transp.",
        **line_kwargs,
    )
    ax.scatter(
        u_vals,
        inv_sim_p,
        marker="o",
        s=LEFT_SIM_MARKER_SIZE**2,
        zorder=10,
        facecolors=color_sim_p,
        edgecolors=color_sim_p,
        linewidths=MARKER_EDGEWIDTH,
        label="Sim. (Post-transp.)",
    )
    if has_m0_sim:
        ax.scatter(
            u_vals,
            inv_sim_p_m0,
            marker="^",
            s=LEFT_SIM_MARKER_SIZE**2,
            zorder=10,
            facecolors=color_m0_marker,
            edgecolors=color_m0_marker,
            linewidths=MARKER_EDGEWIDTH,
            label="Sim. (pre-transp.)",
        )
    style_log_xaxis(ax, u_vals)
    ax.set_xlabel("Transposition rate ($u$)", fontsize=X_AXIS_LABELSIZE, fontweight="bold")
    ax.set_ylabel(r"Dispersion ratio ($\rho$)", fontsize=Y_AXIS_LABELSIZE, fontweight="bold", rotation=90, labelpad=14)
    ax.tick_params(axis="y", labelsize=Y_TICK_LABELSIZE)
    for label in ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.grid(
        SHOW_HORIZONTAL_GRID,
        which="major",
        axis="y",
        linestyle=HORIZONTAL_GRID_LINESTYLE,
        linewidth=HORIZONTAL_GRID_LINEWIDTH,
        alpha=HORIZONTAL_GRID_ALPHA,
    )
    ax.grid(
        SHOW_VERTICAL_GRID,
        which="major",
        axis="x",
        linestyle=VERTICAL_GRID_LINESTYLE,
        linewidth=VERTICAL_GRID_LINEWIDTH,
        alpha=VERTICAL_GRID_ALPHA,
    )
    ax.legend(loc="best", prop={"size": LEFT_LEGEND_FONTSIZE, "weight": "bold"})


# ----- Combined figure: Plot 2 and Plot 3 side by side
inv_app_p = np.array([1.0 / max(p, 1e-12) for p in app_p])
inv_nb_p = np.array([1.0 / max(p, 1e-12) for p in nb_p])
inv_sim_p = np.array([1.0 / max(p, 1e-12) for p in sim_p])
if has_m0_sim:
    inv_sim_p_m0 = np.array([1.0 / max(p, 1e-12) for p in sim_p_m0])

color_red = "#d62728"  # Left panel: stochastic simulations (markers), also used for the pre-transp. triangle marker
color_roze = "#0072B2"  # Roze (2023); colorblind-safe blue, dotted line
color_roze3u = "#2ca02c"  # Roze (2023) 1+3u (post-transposition); green
color_m0 = "#E69F00"  # our own pre-transposition prediction (line); colorblind-safe amber, distinct from the dotted blue Roze line
color_m0_marker = color_red  # pre-transp. simulation triangle marker
color_black = "#000000"  # Left panel: approximate & parametric theory lines
color_sim_p = color_red
color_app = color_black
color_nb = color_black

fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(16, 6), dpi=120)

x_smooth_app, y_smooth_app = smooth_log_spline(x_line_u, inv_app_p[idx_sort_u])
x_smooth_nb, y_smooth_nb = smooth_log_spline(x_line_u, inv_nb_p[idx_sort_u])
x_smooth_roze, y_smooth_roze = smooth_log_spline(x_line_u, rho_roze[idx_sort_u])
x_smooth_roze3u, y_smooth_roze3u = smooth_log_spline(x_line_u, rho_roze_3u[idx_sort_u])
x_smooth_m0, y_smooth_m0 = smooth_log_spline(x_line_u, rho_m0_vals[idx_sort_u])

line_kwargs = dict(
    linewidth=MAIN_LINEWIDTH,
    zorder=3,
    antialiased=True,
    solid_joinstyle="round",
    solid_capstyle="round",
)

plot_left_panel(ax_left)

if has_skew:
    color_skew = "#1f77b4"
    color_exkurt = "#ff7f0e"
    x_smooth_skew, y_smooth_skew = smooth_log_spline(x_line_u, nb_skew[idx_sort_u])
    x_smooth_exkurt, y_smooth_exkurt = smooth_log_spline(x_line_u, nb_exkurt[idx_sort_u])

    ax_right.scatter(
        u_vals,
        sim_skew,
        color=color_skew,
        marker="o",
        s=RIGHT_SIM_MARKER_SIZE**2,
        zorder=6,
        label="Skewness (sim)",
    )
    ax_right.scatter(
        u_vals,
        sim_exkurt,
        color=color_exkurt,
        marker="o",
        s=RIGHT_SIM_MARKER_SIZE**2,
        zorder=6,
        label="Ex. kurtosis (sim)",
    )
    ax_right.plot(
        x_smooth_skew,
        y_smooth_skew,
        color=color_skew,
        linestyle="-",
        linewidth=SECONDARY_LINEWIDTH,
        zorder=3,
        label="Skewness (parametric)",
    )
    ax_right.plot(
        x_smooth_exkurt,
        y_smooth_exkurt,
        color=color_exkurt,
        linestyle="-",
        linewidth=SECONDARY_LINEWIDTH,
        zorder=3,
        label="Ex. kurtosis (parametric)",
    )

    style_log_xaxis(ax_right, u_vals)
    ax_right.set_xlabel("transposition rate ($u$)", fontsize=X_AXIS_LABELSIZE, fontweight="bold")
    ax_right.set_ylabel("Values", fontsize=Y_AXIS_LABELSIZE, fontweight="bold")
    ax_right.tick_params(axis="y", labelsize=Y_TICK_LABELSIZE)
    for label in ax_right.get_yticklabels():
        label.set_fontweight("bold")
    ax_right.grid(
        SHOW_HORIZONTAL_GRID,
        which="major",
        axis="y",
        linestyle=HORIZONTAL_GRID_LINESTYLE,
        linewidth=HORIZONTAL_GRID_LINEWIDTH,
        alpha=HORIZONTAL_GRID_ALPHA,
    )
    ax_right.grid(
        SHOW_VERTICAL_GRID,
        which="major",
        axis="x",
        linestyle=VERTICAL_GRID_LINESTYLE,
        linewidth=VERTICAL_GRID_LINEWIDTH,
        alpha=VERTICAL_GRID_ALPHA,
    )
    ax_right.legend(loc="best", prop={"size": LEGEND_FONTSIZE, "weight": "bold"})
else:
    ax_right.text(
        0.5,
        0.5,
        "Skewness/excess kurtosis columns not found.",
        ha="center",
        va="center",
        transform=ax_right.transAxes,
        fontsize=12,
    )
    ax_right.set_axis_off()

plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
print(f"Figure saved to: {OUTPUT_PATH}")
plt.show()
