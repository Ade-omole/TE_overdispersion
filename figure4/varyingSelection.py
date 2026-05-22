"""Plot dispersion ratio and higher moments versus selection coefficient s."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from scipy import integrate
from scipy.interpolate import make_interp_spline

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"

# Plot styling
LEFT_SIM_MARKER_SIZE = 8
RIGHT_SIM_MARKER_SIZE = 11
MAIN_LINEWIDTH = 5
SECONDARY_LINEWIDTH = 5.0
MARKER_EDGEWIDTH = 3.5

X_TICK_LABELSIZE = 15
Y_TICK_LABELSIZE = 15
LEGEND_FONTSIZE = 16
X_AXIS_LABELSIZE = 18
Y_AXIS_LABELSIZE = 20

SHOW_VERTICAL_GRID = True
VERTICAL_GRID_ALPHA = 0.18
VERTICAL_GRID_LINESTYLE = "-"
VERTICAL_GRID_LINEWIDTH = 0.8

SHOW_HORIZONTAL_GRID = True
HORIZONTAL_GRID_ALPHA = 0.18
HORIZONTAL_GRID_LINESTYLE = "-"
HORIZONTAL_GRID_LINEWIDTH = 0.8

LEFT_LEGEND_X = 1.0
LEFT_LEGEND_Y = 0.09
LEFT_LEGEND_LOC = "lower right"
EVENT_MODEL = "binomial"  # must match varyingSelection.jl
FITNESS_TAG = "exp"

R_MAP = 10.0
U = 1e-2
V = 1e-4
ALPHA_ROZE = 0.0

CSV_SIM = CSV_DIR / f"varyingSelection_{EVENT_MODEL}_{FITNESS_TAG}_u{U}_v{V}_4.csv"

def compute_E1(R: float, u: float, eps: float = 1e-12) -> float:
    """ε₁ for linear map of length R Morgans (Roze appendix)."""
    R = max(R, eps)
    rho_R = u / R
    if R <= 1.0:
        half_plus_rho = 0.5 + rho_R
        if rho_R <= 0.0 or half_plus_rho <= 0.0:
            return 0.0
        return (2.0 / R) * (
            (1.0 + 2.0 * rho_R) * (np.log(half_plus_rho) - np.log(max(rho_R, eps))) - 1.0
        )

    def integrand(x: float) -> float:
        denom = (1.0 - np.exp(-2.0 * x)) / 2.0 + 2.0 * u
        return (R - x) / max(np.finfo(float).eps, denom)

    res, _ = integrate.quad(integrand, 0.0, R, limit=200)
    return (2.0 / (R**2)) * res


def roze_rho_eq14(E1: float, u: float, v: float, alpha: float = 0.0, eps: float = 1e-10) -> float:
    """Roze (2023) eq. (14)."""
    denom = 1.0 - u * E1
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    return 1.0 + (E1 / denom) * ((u + v + alpha) / 2.0)


def roze_rho_eq5_13(E1: float, u: float, v: float, alpha: float = 0.0) -> float:
    """ρ = Var/Mean from n̄ (1 + ε₁ (u+v+α)/2); Roze (2023) eqs. (5) and (13)."""
    return 1.0 + E1 * ((u + v + alpha) / 2.0)


def roze_mean_eq11(s: float, u: float, v: float, alpha: float = 0.0) -> float:
    """Roze (2023) eq.~(11): n̄ ≈ (u − v − α) / β with β = 2s (same β as w(n)=exp(−s n²) in this script)."""
    beta = 2.0 * max(float(s), np.finfo(float).eps)
    num = u - v - alpha
    if num <= 0.0:
        return np.nan
    return num / beta


def roze_variance_from_rho513(mu_roze: float, rho_513: float) -> float:
    """σ² = ρ · n̄ pairing eqs.~(5)+(13) with eq.~(11) mean (ρ from ``roze_rho_eq5_13``)."""
    if not (np.isfinite(mu_roze) and np.isfinite(rho_513)):
        return np.nan
    return rho_513 * mu_roze


def approx_mean_variance(s: float, u: float, v: float):
    """Closed-form approximate equilibrium mean and variance (varyingRsim.jl)."""
    denom_app = 2.0 * s * (2.0 * u + 2.0 * v + 1.0)
    if denom_app <= 0:
        return 0.0, np.nan
    mu_app = max(0.0, ((u - v) * (1.0 - 2.0 * (u - v))) / denom_app)
    sig2_app = max(np.finfo(float).eps, (u - v) / (2.0 * s))
    return mu_app, sig2_app


def dln_omega_dmu(mu: float, s: float, fit: str) -> float:
    return -2.0 * s * mu if fit == "exp" else ((-2.0 * s * mu) / max(np.finfo(float).eps, 1.0 - s * mu**2))


def d2ln_omega_dmu2(mu: float, s: float, fit: str) -> float:
    return -2.0 * s if fit == "exp" else (
        (-2.0 * s - 2.0 * s**2 * mu**2) / max(np.finfo(float).eps, (1.0 - s * mu**2) ** 2)
    )


def nb_equilibrium_moments(
    s: float,
    u: float,
    v: float,
    fit: str = "exp",
    dt: float = 1.0,
    max_gen: int = 80_000,
    tol: float = 1e-7,
    check_every: int = 2_000,
):
    """
    NB parametric closure: same ODE as varyingRsim.jl / varyingSelection.jl.
    Returns (rho, skew, excess_kurtosis, mu, sigma^2) at equilibrium, with
      rho = sigma^2/mu, skew = rho_shape/sqrt(sigma^2), ex. kurt = alpha_shape/sigma^2
    (same as Julia nb_skew_last, nb_exkurt_last after the loop; mu and sigma^2 are the closure state).
    """
    eps = np.finfo(float).eps
    mu_nb = 10.0
    sig2_nb = 10.0
    prev_mu, prev_s2 = mu_nb, sig2_nb

    for gen in range(1, max_gen + 1):
        b1 = dln_omega_dmu(mu_nb, s, fit)
        b2 = d2ln_omega_dmu2(mu_nb, s, fit)
        p_nb = float(np.clip(mu_nb / max(eps, sig2_nb), eps, 1.0 - eps))
        rho_shape = (2.0 - p_nb) / p_nb  # (2-p)/p in closure, not Var/mean
        alpha_nb = (6.0 * (1.0 - p_nb) + p_nb**2) / max(eps, p_nb**2)
        dmu = (u - v) * mu_nb + sig2_nb * b1 + 0.5 * rho_shape * sig2_nb * b2
        ds2 = (
            (2.0 * u + 0.5) * mu_nb
            + (u - v - 0.5) * sig2_nb
            + 0.5 * (b1 + b2 * mu_nb) * sig2_nb
            + 0.25 * (2.0 * b1 + b2) * rho_shape * sig2_nb
            + 0.25 * alpha_nb * b2 * sig2_nb
        )
        mu_nb = max(0.0, mu_nb + dmu * dt)
        sig2_nb = max(eps, sig2_nb + ds2 * dt)

        if gen % check_every == 0:
            rel = max(
                abs(mu_nb - prev_mu) / max(abs(mu_nb), 1.0),
                abs(sig2_nb - prev_s2) / max(abs(sig2_nb), 1.0),
            )
            if rel < tol:
                break
            prev_mu, prev_s2 = mu_nb, sig2_nb

    if mu_nb <= 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    sig2_nb = max(eps, sig2_nb)
    p_nb = float(np.clip(mu_nb / sig2_nb, eps, 1.0 - eps))
    rho_shape = (2.0 - p_nb) / p_nb
    alpha_nb = (6.0 * (1.0 - p_nb) + p_nb**2) / max(eps, p_nb**2)
    rho_disp = sig2_nb / mu_nb
    skew = rho_shape / np.sqrt(sig2_nb)
    exkurt = alpha_nb / sig2_nb
    return (rho_disp, skew, exkurt, mu_nb, sig2_nb)


# Table formatting (used by main only)


def _fmt_rho(x) -> str:
    """Dispersion ratios ρ: four decimal places."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{float(x):.4f}"


def _fmt_two(x) -> str:
    """Non-ρ numeric cells: two decimal places."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):d}"
    return f"{float(x):.2f}"


def _fmt_s_beta(x) -> str:
    """``s`` and ``β = 2s`` can be ≪ 1; two fixed decimals would read as 0.00, so use compact general format."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{float(x):.6g}"


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
    """Format tick labels for log scale with bold math formatting."""
    exponent = int(np.floor(np.log10(value)))
    coeff = value / (10**exponent)
    coeff_rounded = int(round(coeff))
    if np.isclose(coeff, 1.0):
        return rf"$\mathbf{{10^{{{exponent}}}}}$"
    return rf"$\mathbf{{{coeff_rounded} \times 10^{{{exponent}}}}}$"


def style_log_xaxis(ax, x_values):
    """Apply styled log-scale formatting to x-axis."""
    ax.set_xscale("log")
    candidate_ticks = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
    x_min = np.min(x_values)
    x_max = np.max(x_values)
    visible_ticks = candidate_ticks[(candidate_ticks >= x_min) & (candidate_ticks <= x_max)]
    if len(visible_ticks) == 0:
        visible_ticks = candidate_ticks
    ax.xaxis.set_major_locator(FixedLocator(visible_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(_format_log_tick))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", which="major", labelsize=X_TICK_LABELSIZE, length=5, width=1.0, direction="out")


def _fmt_sim(sim: pd.Series | None, key: str) -> str:
    if sim is None or key not in sim.index:
        return "—"
    v = sim[key]
    if pd.isna(v):
        return "—"
    if key == "rho_sim":
        return _fmt_rho(float(v))
    return _fmt_two(float(v))


def _sim_at_s(df: pd.DataFrame, s: float, rtol: float = 1e-6, atol: float = 1e-12) -> pd.Series | None:
    """First simulation row matching ``s`` (float-safe), else None."""
    ss = df["_s"].to_numpy(dtype=float)
    j = int(np.argmin(np.abs(ss - s)))
    if abs(ss[j] - s) <= atol + rtol * max(abs(s), 1e-30):
        return df.iloc[j]
    return None


def main():
    # Keep in sync with s_list in varyingSelection.jl.
    s_values = np.sort(
        np.unique(
            np.array([3e-5, 1e-4, 3e-4, 5e-4, 1e-3, 2e-3], dtype=float),
        )
    )
    beta_values = 2.0 * s_values

    mu_app = np.empty_like(beta_values)
    sig2_app = np.empty_like(beta_values)
    rho_app = np.empty_like(beta_values)
    rho_nb = np.empty_like(beta_values)
    mu_nb = np.empty_like(beta_values)
    sig2_nb = np.empty_like(beta_values)
    nb_skew = np.empty_like(beta_values)
    nb_exkurt = np.empty_like(beta_values)
    rho_roze_14 = np.empty_like(beta_values)
    rho_roze_513 = np.empty_like(beta_values)
    mu_roze = np.empty_like(beta_values)
    sig2_roze = np.empty_like(beta_values)

    E1 = compute_E1(R_MAP, U)

    for i, s in enumerate(s_values):
        mu_a, var_a = approx_mean_variance(s, U, V)
        mu_app[i] = mu_a
        sig2_app[i] = var_a
        rho_app[i] = var_a / mu_a if mu_a > 0 and np.isfinite(var_a) else np.nan
        rnb, sk, ek, mnb, s2nb = nb_equilibrium_moments(s, U, V, fit="exp")
        rho_nb[i] = rnb
        mu_nb[i] = mnb
        sig2_nb[i] = s2nb
        nb_skew[i] = sk
        nb_exkurt[i] = ek
        rho_roze_14[i] = roze_rho_eq14(E1, U, V, ALPHA_ROZE)
        rho_roze_513[i] = roze_rho_eq5_13(E1, U, V, ALPHA_ROZE)
        mu_roze[i] = roze_mean_eq11(s, U, V, ALPHA_ROZE)
        sig2_roze[i] = roze_variance_from_rho513(mu_roze[i], rho_roze_513[i])

    df_sim = None
    if CSV_SIM.is_file():
        df_sim = pd.read_csv(CSV_SIM)
        df_sim = df_sim.copy()
        df_sim["_s"] = df_sim["s"].astype(float) if "s" in df_sim.columns else (df_sim["beta"].astype(float) / 2.0)
    else:
        print(f"No {CSV_SIM} — simulation columns will be empty (run: julia varyingSelection.jl)")

    print("\nRoze eq.~(11)  n̄ = (u−v−α)/β  and  σ² = ρ₅,₁₃·n̄  (β = 2s; ρ₅,₁₃ from eqs.~(5)+(13)):")
    for i, s in enumerate(s_values):
        print(
            f"  s={float(s):.6g}  β={beta_values[i]:.6g}  μ_Roze={mu_roze[i]:.2f}  σ²_Roze={sig2_roze[i]:.2f}"
        )

    col_labels = [
        "s",
        "β=2s",
        "μ app",
        "σ² app",
        "ρ app",
        "μ NB",
        "σ² NB",
        "ρ NB",
        "μ Roze",
        "σ² Roze",
        "ρ R14",
        "ρ R513",
        "ρ sim",
        "μ sim",
        "σ² sim",
        "skew NB",
        "ex.κ NB",
        "skew sim",
        "ex.κ sim",
    ]
    cell_text: list[list[str]] = []
    for i, s in enumerate(s_values):
        s = float(s)
        sim = _sim_at_s(df_sim, s) if df_sim is not None else None
        cell_text.append(
            [
                _fmt_s_beta(s),
                _fmt_s_beta(beta_values[i]),
                _fmt_two(mu_app[i]),
                _fmt_two(sig2_app[i]),
                _fmt_rho(rho_app[i]),
                _fmt_two(mu_nb[i]),
                _fmt_two(sig2_nb[i]),
                _fmt_rho(rho_nb[i]),
                _fmt_two(mu_roze[i]),
                _fmt_two(sig2_roze[i]),
                _fmt_rho(rho_roze_14[i]),
                _fmt_rho(rho_roze_513[i]),
                _fmt_sim(sim, "rho_sim"),
                _fmt_sim(sim, "sim_μ"),
                _fmt_sim(sim, "sim_σ²"),
                _fmt_two(nb_skew[i]),
                _fmt_two(nb_exkurt[i]),
                _fmt_sim(sim, "sim_skew"),
                _fmt_sim(sim, "sim_exkurt"),
            ]
        )

    nrows = len(s_values) + 1
    fig_h = min(22.0, max(6.0, 0.52 * nrows))
    fig_w = 26.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=120)
    ax.axis("off")
    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="upper center",
        cellLoc="right",
        bbox=[0.02, 0.02, 0.96, 0.92],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.0, 1.65)
    for key, cell in tbl.get_celld().items():
        r, c = key
        if r == 0:
            cell.set_facecolor("#e0e0e0")
            cell.set_text_props(weight="bold")
            cell.set_height(0.055)
        else:
            cell.set_facecolor("#f8f8f8" if r % 2 == 0 else "#ffffff")

    title = (
        rf"$\rho$ and moments by $s$ — $R={R_MAP:g}$, $u={U:g}$, $v={V:g}$, $\alpha={ALPHA_ROZE:g}$, "
        rf"$\varepsilon_1={E1:.5g}$"
    )
    fig.suptitle(title, fontsize=11, y=0.995)

    out_table = FIG_DIR / "varyingSelection_table.pdf"
    fig.savefig(out_table, bbox_inches="tight")
    print(f"Saved {out_table}")

    color_red = "#d62728"
    color_roze = "#9467bd"
    color_black = "#000000"
    color_skew = "#1f77b4"
    color_exkurt = "#ff7f0e"

    line_kwargs = dict(
        linewidth=MAIN_LINEWIDTH,
        zorder=3,
        antialiased=True,
        solid_joinstyle="round",
        solid_capstyle="round",
    )

    idx_sort_s = np.argsort(s_values)
    x_line_s = s_values[idx_sort_s]
    x_smooth_app, y_smooth_app = smooth_log_spline(x_line_s, rho_app[idx_sort_s])
    x_smooth_nb, y_smooth_nb = smooth_log_spline(x_line_s, rho_nb[idx_sort_s])
    x_smooth_roze14, y_smooth_roze14 = smooth_log_spline(x_line_s, rho_roze_14[idx_sort_s])
    x_smooth_skew, y_smooth_skew = smooth_log_spline(x_line_s, nb_skew[idx_sort_s])
    x_smooth_exkurt, y_smooth_exkurt = smooth_log_spline(x_line_s, nb_exkurt[idx_sort_s])

    fig_combined, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
    ax_left.plot(x_smooth_app, y_smooth_app, color=color_black, linestyle="--", label="Approximate", **line_kwargs)
    ax_left.plot(x_smooth_nb, y_smooth_nb, color=color_black, linestyle="-", label="Parametric", **line_kwargs)
    ax_left.plot(x_smooth_roze14, y_smooth_roze14, color=color_roze, linestyle="-", label="Roze (2023)", **line_kwargs)

    if df_sim is not None:
        ax_left.scatter(
            df_sim["_s"],
            df_sim["rho_sim"],
            marker="o",
            s=LEFT_SIM_MARKER_SIZE**2,
            zorder=10,
            facecolors=color_red,
            edgecolors=color_red,
            linewidths=MARKER_EDGEWIDTH,
            label="Simulation",
        )

    style_log_xaxis(ax_left, s_values)
    ax_left.set_xlabel("Selection coefficient ($s$)", fontsize=X_AXIS_LABELSIZE, fontweight="bold")
    ax_left.set_ylabel(r"Dispersion ratio ($\rho$)", fontsize=Y_AXIS_LABELSIZE, fontweight="bold", rotation=90, labelpad=14)
    ax_left.tick_params(axis="y", labelsize=Y_TICK_LABELSIZE)
    for label in ax_left.get_yticklabels():
        label.set_fontweight("bold")
    ax_left.grid(
        SHOW_HORIZONTAL_GRID,
        which="major",
        axis="y",
        linestyle=HORIZONTAL_GRID_LINESTYLE,
        linewidth=HORIZONTAL_GRID_LINEWIDTH,
        alpha=HORIZONTAL_GRID_ALPHA,
    )
    ax_left.grid(
        SHOW_VERTICAL_GRID,
        which="major",
        axis="x",
        linestyle=VERTICAL_GRID_LINESTYLE,
        linewidth=VERTICAL_GRID_LINEWIDTH,
        alpha=VERTICAL_GRID_ALPHA,
    )
    ax_left.legend(
        loc=LEFT_LEGEND_LOC,
        bbox_to_anchor=(LEFT_LEGEND_X, LEFT_LEGEND_Y),
        prop={"size": LEGEND_FONTSIZE, "weight": "bold"},
    )

    has_sim_moments = df_sim is not None and "sim_skew" in df_sim.columns and "sim_exkurt" in df_sim.columns

    if has_sim_moments:
        ax_right.scatter(
            df_sim["_s"],
            df_sim["sim_skew"],
            color=color_skew,
            marker="o",
            s=RIGHT_SIM_MARKER_SIZE**2,
            zorder=6,
            label="Skewness (sim)",
        )
        ax_right.scatter(
            df_sim["_s"],
            df_sim["sim_exkurt"],
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

    style_log_xaxis(ax_right, s_values)
    ax_right.set_xlabel("Selection coefficient ($s$)", fontsize=X_AXIS_LABELSIZE, fontweight="bold")
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

    plt.tight_layout()
    out_combined = FIG_DIR / "varyingSelection.pdf"
    fig_combined.savefig(out_combined, bbox_inches="tight")
    print(f"Saved {out_combined}")

    print(f"ε₁(R={R_MAP}, u={U}) = {E1:.4f}  (Roze ρ from eq.~(14) and (5)+(13) independent of β here)")
    _figures = [fig, fig_combined]
    if os.environ.get("MPLBACKEND", "").lower() != "agg":
        plt.show()
    else:
        for f in _figures:
            plt.close(f)


if __name__ == "__main__":
    main()
