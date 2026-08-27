"""Plot dispersion ratio under varying fitness curvature (Appendix Fig. E.1)."""

from pathlib import Path

import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"

# Fitness function: "exp" for w(n) = exp(-s n^gamma)
FITNESS_TYPE = "exp"  # "exp" or "synergistic"

# Parameters, matching Fig. 3
N_VAL = 1e4
S_VAL = 1.0 / N_VAL  # s = 1/N
V_VAL = 1e-4
U_MIN, U_MAX, N_U = 5e-4, 0.032, 25  # matches varyingFitnesssim.jl's u_values range
GAMMAS = [1.5, 2, 2.5, 3, 4, 7, 10]

# Style, following varyingRateFig.py
MAIN_LINEWIDTH = 5  # theory curves in all three panels
SIM_MARKER_SIZE = 7  # simulation markers (scatter takes the squared size)
TICK_LABELSIZE = 15  # x- and y-axis tick labels
AXIS_LABELSIZE = 18  # axis labels
TITLE_FONTSIZE = 18  # panel titles
LEGEND_FONTSIZE = 13  # legend entries, its title, and the "Simulation" heading
END_LABEL_FONTSIZE = 12  # direct end-of-line gamma labels

def beta1_exp(mu, s, gamma):
    """beta1 = d/dmu ln w(mu) for w(n) = exp(-s n^gamma)."""
    return -s * gamma * mu ** (gamma - 1.0)


def beta2_exp(mu, s, gamma):
    """beta2 = d^2/dmu^2 ln w(mu) for w(n) = exp(-s n^gamma)."""
    return -s * gamma * (gamma - 1.0) * mu ** (gamma - 2.0)


def domain_margin_exp(mu, s, gamma):
    """w(n)=exp(-s n^gamma) is defined for all mu>0: no domain boundary."""
    return None


def w_exp(n, s, gamma):
    return np.exp(-s * n**gamma)


def beta1_synergistic(mu, s, gamma):
    """beta1 = d/dmu ln w(mu) for w(n) = 1 - s n^gamma."""
    f = 1.0 - s * mu**gamma
    f_prime = -s * gamma * mu ** (gamma - 1.0)
    return f_prime / f


def beta2_synergistic(mu, s, gamma):
    """beta2 = d^2/dmu^2 ln w(mu) for w(n) = 1 - s n^gamma."""
    f = 1.0 - s * mu**gamma
    f_prime = -s * gamma * mu ** (gamma - 1.0)
    f_double_prime = -s * gamma * (gamma - 1.0) * mu ** (gamma - 2.0)
    b1 = f_prime / f
    return f_double_prime / f - b1**2


def domain_margin_synergistic(mu, s, gamma):
    """w(n) = 1 - s n^gamma is only defined (and positive) while s*mu^gamma < 1;"""
    return 1.0 - s * mu**gamma


def w_synergistic(n, s, gamma):
    return 1.0 - s * n**gamma


def n_pole_synergistic(s, gamma):
    """n at which w(n) = 1 - s n^gamma hits zero."""
    return (1.0 / s) ** (1.0 / gamma)


FITNESS_MODELS = {
    "exp": dict(
        beta1=beta1_exp,
        beta2=beta2_exp,
        domain_margin=domain_margin_exp,
        w=w_exp,
        n_pole=None,  # no finite domain boundary
        legend_title=r"$w(n)=\exp(-sn^{\gamma})$",
        tag="exp",
    ),
    "synergistic": dict(
        beta1=beta1_synergistic,
        beta2=beta2_synergistic,
        domain_margin=domain_margin_synergistic,
        w=w_synergistic,
        n_pole=n_pole_synergistic,
        legend_title=r"$w(n)=1-sn^{\gamma}$",
        tag="synergistic",
    ),
}
_model = FITNESS_MODELS[FITNESS_TYPE]
beta1, beta2, domain_margin = _model["beta1"], _model["beta2"], _model["domain_margin"]
fitness_w, n_pole = _model["w"], _model["n_pole"]
OUTPUT_PATH = FIG_DIR / f"varyingFitness_{_model['tag']}.pdf"

# Simulation overlay (from varyingFitnesssim.jl)
SIM_TAG = "exp"
PRE_SIM_GAMMAS = [2.5, 4] if FITNESS_TYPE == "exp" else []
POST_SIM_GAMMAS = [2.5, 4] if FITNESS_TYPE == "exp" else []


def nb_shape_terms(mu, sigma2):
    """vartheta and alpha (Eqs. eq:nb_skew, eq:nb_exkurt), from p = mu/sigma2."""
    p = mu / sigma2
    vartheta = (2.0 - p) / p
    alpha = (6.0 * (1.0 - p) + p**2) / p**2
    return vartheta, alpha


def nb_closure_rhs_log(t, log_state, u, v, s, gamma):
    """Eq. ODE_Nb, parameterized in log(mu), log(sigma2) so mu, sigma2 stay"""
    mu, sigma2 = np.exp(log_state)
    b1, b2 = beta1(mu, s, gamma), beta2(mu, s, gamma)
    vartheta, alpha = nb_shape_terms(mu, sigma2)

    dmu = (u - v) * mu + sigma2 * b1 + 0.5 * vartheta * sigma2 * b2
    dsigma2 = (
        (2.0 * u + 0.5) * mu
        + (u - v - 0.5) * sigma2
        + 0.5 * (b1 + b2 * mu) * sigma2
        + 0.25 * (2.0 * b1 + b2) * vartheta * sigma2
        + 0.25 * alpha * b2 * sigma2
    )
    return [dmu / mu, dsigma2 / sigma2]


T_FINAL = 30000.0  # upper bound on integration time, in case steady state is never reached
DOMAIN_MARGIN_TOL = 0.05  # stop the integration once w(n) gets this close to its domain boundary
CONVERGED_RATE_TOL = 1e-9  # stop once relative rates of change drop below this (steady state reached)


def _domain_boundary_event(t, log_state, u, v, s, gamma):
    mu = np.exp(log_state[0])
    return domain_margin(mu, s, gamma) - DOMAIN_MARGIN_TOL


_domain_boundary_event.terminal = True
_domain_boundary_event.direction = -1


def _converged_event(t, log_state, u, v, s, gamma):
    """Fires once the relative rates of change (d log mu/dt, d log sigma2/dt) both"""
    d_log_mu, d_log_sigma2 = nb_closure_rhs_log(t, log_state, u, v, s, gamma)
    return max(abs(d_log_mu), abs(d_log_sigma2)) - CONVERGED_RATE_TOL


_converged_event.terminal = True
_converged_event.direction = -1


class NoValidEquilibrium(Exception):
    """Raised when no valid, overdispersed NB-closure equilibrium exists at this"""


def solve_nb_equilibrium(u, v, s, gamma, mu0, sigma2_0):
    """NB closure equilibrium (mu_nb, sigma2_nb), found by integrating Eq. ODE_Nb"""
    log_state0 = [np.log(mu0), np.log(sigma2_0)]
    has_domain_boundary = domain_margin(mu0, s, gamma) is not None
    events = [_converged_event, _domain_boundary_event] if has_domain_boundary else [_converged_event]
    sol = solve_ivp(
        nb_closure_rhs_log,
        [0.0, T_FINAL],
        log_state0,
        args=(u, v, s, gamma),
        method="Radau",
        rtol=1e-10,
        atol=1e-13,
        max_step=200.0,  # avoids the solver stalling/hanging near a domain boundary
        events=events,
    )
    if not sol.success:
        raise RuntimeError(f"NB equilibrium integration failed for u={u}, gamma={gamma}: {sol.message}")
    domain_event_idx = 1 if has_domain_boundary else None
    if domain_event_idx is not None and len(sol.t_events[domain_event_idx]) > 0:
        raise NoValidEquilibrium(f"w(n) approached its domain boundary for u={u}, gamma={gamma}")
    mu, sigma2 = np.exp(sol.y[:, -1])
    if sigma2 <= mu or mu <= 0:
        raise NoValidEquilibrium(f"NB equilibrium not overdispersed for u={u}, gamma={gamma}: mu={mu}, sigma2={sigma2}")
    return mu, sigma2


def pretransposition_moments(mu, sigma2, s, gamma):
    """Pre-transposition/excision E[m'] and Var[m'] (Eqs. eq:premean_nb, eq:prevar_nb)."""
    b1, b2 = beta1(mu, s, gamma), beta2(mu, s, gamma)
    vartheta, alpha = nb_shape_terms(mu, sigma2)

    e_m = mu + b1 * sigma2 + 0.5 * vartheta * b2 * sigma2
    var_m = (
        0.5 * (mu + sigma2)
        + 0.5 * b1 * (1.0 + vartheta) * sigma2
        + 0.25 * b2 * ((vartheta + alpha) * sigma2 + 2.0 * sigma2**2)
    )
    return e_m, var_m


def sweep_gamma(gamma, u_vals, v, s):
    """Solve the NB equilibrium across u_vals for one gamma, using continuation"""
    # Start from a rare invader: mu0 = sigma2_0 = 0.5
    mu_guess, sigma2_guess = 0.5, 0.5

    rho_post = np.full_like(u_vals, np.nan)
    rho_pre = np.full_like(u_vals, np.nan)
    for i, u in enumerate(u_vals):
        try:
            mu, sigma2 = solve_nb_equilibrium(u, v, s, gamma, mu_guess, sigma2_guess)
        except NoValidEquilibrium as exc:
            # No overdispersed equilibrium here: leave the rest as NaN
            print(f"  [gamma={gamma:g}] stopping this curve at u={u:.5f}: {exc}")
            break

        rho_post[i] = sigma2 / mu
        e_m, var_m = pretransposition_moments(mu, sigma2, s, gamma)
        rho_pre[i] = var_m / e_m

        mu_guess, sigma2_guess = mu, sigma2  # continuation seed for next u
    return rho_pre, rho_post


# Sweep u for every gamma
u_vals = np.geomspace(U_MIN, U_MAX, N_U)
results = {gamma: sweep_gamma(gamma, u_vals, V_VAL, S_VAL) for gamma in GAMMAS}


def load_sim_data():
    """Concatenate all varyingFitness_sim_<SIM_TAG>_gamma*.csv files (one per"""
    pattern = str(CSV_DIR / f"varyingFitness_sim_{SIM_TAG}_gamma*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


sim_df = load_sim_data()
if sim_df is None:
    print(f"No simulation CSVs found matching varyingFitness_sim_{SIM_TAG}_gamma*.csv -- theory-only plot.")


# Plotting
def smooth_log_spline(x_sorted, y_sorted, n_interp=400):
    """Cubic spline in log10(x), dropping any trailing NaNs (a fitness model"""
    finite = np.isfinite(y_sorted)
    x_sorted, y_sorted = x_sorted[finite], y_sorted[finite]
    log_x = np.log10(x_sorted)
    if len(log_x) < 4:
        return x_sorted, y_sorted
    spline = make_interp_spline(log_x, y_sorted, k=3)
    log_x_fine = np.linspace(log_x[0], log_x[-1], n_interp)
    return 10**log_x_fine, spline(log_x_fine)


def _format_log_tick(value, _pos):
    exponent = int(np.floor(np.log10(value)))
    coeff = value / (10**exponent)
    if np.isclose(coeff, 1.0):
        return rf"$\mathbf{{10^{{{exponent}}}}}$"
    return rf"$\mathbf{{{coeff:.0f} \times 10^{{{exponent}}}}}$"


def style_log_xaxis(ax, x_min, x_max, candidate_ticks=None):
    """Log-scale x-axis with bold power-of-ten tick labels. If candidate_ticks"""
    ax.set_xscale("log")
    if candidate_ticks is None:
        lo, hi = int(np.floor(np.log10(x_min))), int(np.ceil(np.log10(x_max)))
        candidate_ticks = 10.0 ** np.arange(lo, hi + 1)
    candidate_ticks = np.asarray(candidate_ticks)
    visible = candidate_ticks[(candidate_ticks >= x_min) & (candidate_ticks <= x_max)]
    ax.xaxis.set_major_locator(FixedLocator(visible))
    ax.xaxis.set_major_formatter(FuncFormatter(_format_log_tick))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", which="major", labelsize=TICK_LABELSIZE)


def bold_tick_labels(ax, axis="y"):
    """Bold the numeric tick labels (the log x-axis labels above are already"""
    labels = ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()
    for label in labels:
        label.set_fontweight("bold")


# Paul Tol's qualitative palette
GAMMA_COLORS = [
    "#332288",  # indigo
    "#88CCEE",  # cyan
    "#117733",  # green
    "#DDCC77",  # sand
    "#CC6677",  # rose
    "#AA4499",  # purple
    "#882255",  # wine
]

# Marker shape per gamma
SIM_COLOR = "red"
GAMMA_MARKERS = {2.5: "^", 4: "o"}


def add_end_labels(ax, entries, min_frac=0.045):
    """Direct end-of-line labels, colored to match each curve, staggered"""
    ax.figure.canvas.draw()  # ensure ylim reflects the plotted data
    y_lo, y_hi = ax.get_ylim()
    min_gap = min_frac * (y_hi - y_lo)

    # Push overlapping labels apart, outward from the median
    ordered = sorted(range(len(entries)), key=lambda i: entries[i][0])
    y_placed = [entries[i][0] for i in ordered]
    mid = len(y_placed) // 2
    for i in range(mid + 1, len(y_placed)):
        if y_placed[i] - y_placed[i - 1] < min_gap:
            y_placed[i] = y_placed[i - 1] + min_gap
    for i in range(mid - 1, -1, -1):
        if y_placed[i + 1] - y_placed[i] < min_gap:
            y_placed[i] = y_placed[i + 1] - min_gap

    for idx, y_label in zip(ordered, y_placed):
        y_end, x_end, gamma, color = entries[idx]
        ax.annotate(
            rf"$\gamma={gamma:g}$",
            xy=(x_end, y_end),
            xytext=(x_end, y_label),
            textcoords="data",
            xycoords="data",
            color=color,
            fontsize=END_LABEL_FONTSIZE,
            fontweight="bold",
            va="center",
            annotation_clip=False,
            arrowprops=dict(arrowstyle="-", color=color, alpha=0.5, linewidth=0.8, shrinkA=0, shrinkB=2),
        )


# w(n) itself, for the leftmost panel
FITNESS_DISPLAY_N0 = 25.0
N_MIN_PLOT = 0.0
N_MAX_PLOT = FITNESS_DISPLAY_N0 if n_pole is not None else 1.3 * FITNESS_DISPLAY_N0


def fitness_curve(gamma):
    """n and w(n) arrays for one gamma, using s_gamma=FITNESS_DISPLAY_N0^-gamma"""
    s_display = FITNESS_DISPLAY_N0 ** (-gamma)
    n = np.linspace(N_MIN_PLOT, N_MAX_PLOT, 400)
    return n, np.clip(fitness_w(n, s_display, gamma), 0.0, None)


fig, (ax_fit, ax_pre, ax_post) = plt.subplots(1, 3, figsize=(21, 6), dpi=120)

pre_end_entries, post_end_entries = [], []
for gamma, color in zip(GAMMAS, GAMMA_COLORS):
    label = rf"$\gamma = {gamma:g}$"

    n_fit, w_fit = fitness_curve(gamma)
    ax_fit.plot(n_fit, w_fit, color=color, linewidth=MAIN_LINEWIDTH, label=label)

    rho_pre, rho_post = results[gamma]
    x_pre, pre_s = smooth_log_spline(u_vals, rho_pre)
    x_post, post_s = smooth_log_spline(u_vals, rho_post)
    ax_pre.plot(x_pre, pre_s, color=color, linewidth=MAIN_LINEWIDTH, label=label)
    ax_post.plot(x_post, post_s, color=color, linewidth=MAIN_LINEWIDTH, label=label)
    pre_end_entries.append((pre_s[-1], x_pre[-1], gamma, color))
    post_end_entries.append((post_s[-1], x_post[-1], gamma, color))

    if sim_df is not None and gamma in GAMMA_MARKERS:
        sim_g = sim_df[sim_df["gamma"] == gamma]
        scatter_kwargs = dict(
            marker=GAMMA_MARKERS[gamma], s=SIM_MARKER_SIZE**2, facecolors=SIM_COLOR,
            edgecolors="none", zorder=6,
        )
        # Mean of exactly 0 means stochastic extinction
        if gamma in PRE_SIM_GAMMAS:
            sim_pre = sim_g[sim_g["sim_mean_pre"] > 0]
            if not sim_pre.empty:
                ax_pre.scatter(sim_pre["u"], sim_pre["sim_rho_pre"], **scatter_kwargs)
        if gamma in POST_SIM_GAMMAS:
            sim_post = sim_g[sim_g["sim_mean_post"] > 0]
            if not sim_post.empty:
                ax_post.scatter(sim_post["u"], sim_post["sim_rho_post"], **scatter_kwargs)

ax_fit.set_xlim(N_MIN_PLOT, N_MAX_PLOT)
ax_fit.set_xlabel("Copy number ($n$)", fontsize=AXIS_LABELSIZE, fontweight="bold")
ax_fit.set_ylabel(r"Fitness ($w(n)$)", fontsize=AXIS_LABELSIZE, fontweight="bold")
ax_fit.set_title("Fitness function", fontsize=TITLE_FONTSIZE, fontweight="bold")
ax_fit.set_ylim(-0.03, 1.05)
ax_fit.tick_params(axis="both", which="major", labelsize=TICK_LABELSIZE)
ax_fit.grid(True, which="major", linestyle="-", linewidth=0.8, alpha=0.18)
bold_tick_labels(ax_fit, "x")
bold_tick_labels(ax_fit, "y")

for ax, title in ((ax_pre, "Pre-transposition/excision"), (ax_post, "Post-transposition/excision")):
    style_log_xaxis(ax, U_MIN, U_MAX, candidate_ticks=[1e-3, 3e-3, 1e-2, 3e-2])
    ax.set_xlabel("Transposition rate ($u$)", fontsize=AXIS_LABELSIZE, fontweight="bold")
    ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
    ax.tick_params(axis="y", labelsize=TICK_LABELSIZE)
    ax.grid(True, which="major", linestyle="-", linewidth=0.8, alpha=0.18)
    ax.margins(x=0.12)  # room for the end-of-line gamma labels
    bold_tick_labels(ax, "y")

ax_pre.set_ylabel(r"Dispersion ratio ($\rho$)", fontsize=AXIS_LABELSIZE, fontweight="bold")

# Vertical offset of the "Simulation" legend heading
SIM_LABEL_Y_OFFSET = 0.85

# Legend border padding, in font-size units
LEGEND_BORDER_PAD = 0.4

# Single legend box: theory lines, then the simulation markers
legend_handles, legend_labels = ax_pre.get_legend_handles_labels()
if sim_df is not None:
    # Placeholder reuses a real gamma label to fix the row width
    placeholder_label = legend_labels[-1]
    blank_handle = Line2D([0], [0], color="none", label=placeholder_label)
    sim_handles = [
        Line2D(
            [0], [0], marker=GAMMA_MARKERS[g], color="none", markerfacecolor=SIM_COLOR,
            markeredgecolor="none", markersize=SIM_MARKER_SIZE, label=rf"$\gamma={g:g}$",
        )
        for g in PRE_SIM_GAMMAS
    ]
    legend_handles += [blank_handle] + sim_handles
    legend_labels += [placeholder_label] + [rf"$\gamma={g:g}$" for g in PRE_SIM_GAMMAS]
    # Located by position, not by text
    placeholder_row = len(legend_labels) - len(PRE_SIM_GAMMAS) - 1

pre_legend = ax_pre.legend(
    legend_handles,
    legend_labels,
    loc="upper left",
    prop={"size": LEGEND_FONTSIZE, "weight": "bold"},
    title="Theory",  # the fitness function itself is given in the figure caption
    title_fontproperties={"size": LEGEND_FONTSIZE, "weight": "bold"},
    handlelength=1.2,
    handletextpad=0.5,
)
pre_legend.get_title().set_ha("left")
pre_legend._legend_box.align = "left"

if sim_df is not None:
    fig.canvas.draw()  # finalize layout so window extents below are accurate
    title_x0 = pre_legend.get_title().get_window_extent().transformed(ax_pre.transAxes.inverted()).x0

    sim_text = pre_legend.get_texts()[placeholder_row]
    sim_bb = sim_text.get_window_extent().transformed(ax_pre.transAxes.inverted())
    # Place "Simulation" above the centre of its blank row
    y_center = (sim_bb.y0 + sim_bb.y1) / 2.0 + SIM_LABEL_Y_OFFSET * sim_bb.height
    sim_text.set_alpha(0)

    sim_label = ax_pre.text(
        title_x0, y_center, "Simulation", transform=ax_pre.transAxes,
        fontsize=LEGEND_FONTSIZE, fontweight="bold", ha="left", va="center", zorder=10,
    )

    # Painted on top of the legend, not laid out by it
    fig.canvas.draw()
    spacer = ax_pre.text(0, 0, " ", fontsize=LEGEND_FONTSIZE, fontweight="bold")
    space_w = spacer.get_window_extent().width
    spacer.remove()
    inner_right = (
        pre_legend.get_window_extent().x1
        - LEGEND_BORDER_PAD * LEGEND_FONTSIZE * fig.dpi / 72.0
    )
    overhang = sim_label.get_window_extent().x1 - inner_right
    if overhang > 0 and space_w > 0:
        sim_text.set_text(" " * int(np.ceil(overhang / space_w)) + placeholder_label)
        fig.canvas.draw()

add_end_labels(ax_pre, pre_end_entries)
add_end_labels(ax_post, post_end_entries)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
print(f"Figure saved to: {OUTPUT_PATH}")
plt.show()
