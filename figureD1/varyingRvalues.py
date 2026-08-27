"""Print summary tables for the varying-R sweep (simulation vs theory)."""

from pathlib import Path

import pandas as pd
import numpy as np
from scipy import integrate

FIG_DIR = Path(__file__).resolve().parent
TE_ROOT = FIG_DIR.parent
CSV_DIR = TE_ROOT / "csv_files"

# ----- CSV configuration (edit here)
CSV_BASENAME = "sweep_R_topright_burnin_binomial_pre"  # used in the main text for supp value table (revision2)

FITNESS = "exp"                 # "exp" or "syn"

# ----- Life-cycle stage for our approximate and NB closure rows
LIFE_STAGE = "pre"             # "post" or "pre"

csv_path = CSV_DIR / f"{CSV_BASENAME}_{FITNESS}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run varyingRsim.jl first, or change CSV_BASENAME / FITNESS in the CSV configuration section.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

R = df["R"].values

# Columns: sim_μ, app_μ, nb_μ, sim_σ², app_σ², nb_σ², sim_p, app_p, nb_p
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

# Pre-transposition/excision (m') simulation columns, written by varyingRsim.jl
has_sim_pre = "sim_μ_m0" in df.columns
if has_sim_pre:
    sim_mean_m0 = df["sim_μ_m0"].values
    sim_var_m0 = df["sim_σ²_m0"].values
    sim_skew_m0 = df["sim_skew_m0"].values if "sim_skew_m0" in df.columns else None
    sim_exkurt_m0 = df["sim_exkurt_m0"].values if "sim_exkurt_m0" in df.columns else None

# ----- Apply the life-cycle stage to the approximate and NB closure rows
u_arr = df["u"].values
v_arr = df["v"].values
# s is not in the CSV; recover it from sigma^2_app = (u-v)/(2s)
s_arr = (u_arr - v_arr) / (2.0 * app_var)

if LIFE_STAGE == "post":
    app_mean_st, app_var_st = app_mean, app_var
    nb_mean_st, nb_var_st = nb_mean, nb_var
    nb_skew_st = nb_skew if has_skew else None
    nb_exkurt_st = nb_exkurt if has_skew else None
    show_nb_shape = has_skew

    sim_mean_st, sim_var_st = sim_mean, sim_var
    sim_skew_st = sim_skew if has_skew else None
    sim_exkurt_st = sim_exkurt if has_skew else None
    show_sim_shape = has_skew

elif LIFE_STAGE == "pre":
    net = u_arr - v_arr

    # Approximate closure, Suppl. S4.2.
    factor_app = 1.0 - 2.0 * s_arr * app_var
    app_mean_st = app_mean * factor_app
    app_var_st = 0.5 * (app_mean + app_var) * factor_app

    # Parametric (NB) closure, Eq. (23), unexpanded in mu and sigma^2.
    nb_mean_st = nb_mean * (1.0 - net)
    nb_var_st = 0.5 * (
        nb_mean + nb_var * (1.0 - 3.0 * net) + s_arr * nb_var**2 / nb_mean
    )

    # Higher moments at this stage, Eq. (24) of the main text. The CSV
    eps = 3.0 * u_arr + v_arr
    nb_skew_st = (1.0 + 3.0 * u_arr + v_arr) / np.sqrt(nb_var_st)
    nb_exkurt_st = (1.0 + 3.0 * eps + 1.5 * eps**2) / nb_var_st
    show_nb_shape = True

    # Simulation at the same stage, from varyingRsim.jl's m0 columns.
    if has_sim_pre:
        sim_mean_st, sim_var_st = sim_mean_m0, sim_var_m0
        sim_skew_st, sim_exkurt_st = sim_skew_m0, sim_exkurt_m0
        show_sim_shape = sim_skew_m0 is not None
    else:
        print("Warning: no sim_μ_m0 columns in this CSV; the simulation row is "
              "still post-transposition/excision. Re-run varyingRsim.jl.")
        sim_mean_st, sim_var_st = sim_mean, sim_var
        sim_skew_st = sim_skew if has_skew else None
        sim_exkurt_st = sim_exkurt if has_skew else None
        show_sim_shape = has_skew

else:
    raise ValueError(f"LIFE_STAGE must be 'post' or 'pre', got {LIFE_STAGE!r}")

# Roze (2023): mean eq. (11), variance eqs. (5)+(13)
ALPHA_ROZE = 0.0
U_ROZE = 0.01
V_ROZE = U_ROZE / 10.0
BETA_ROZE = U_ROZE / 10.0


def compute_E1(R, u, eps=1e-12):
    """E₁ (ε₁) for linear genetic map of length R Morgans (Roze). R, u scalars."""
    R = max(R, eps)
    rho_R = u / R
    if R <= 1.0:
        half_plus_rho = 0.5 + rho_R
        if rho_R <= 0.0 or half_plus_rho <= 0.0:
            return 0.0
        return (2.0 / R) * ((1.0 + 2.0 * rho_R) * (np.log(half_plus_rho) - np.log(max(rho_R, eps))) - 1.0)

    def integrand(x):
        denom = (1.0 - np.exp(-2.0 * x)) / 2.0 + 2.0 * u
        return (R - x) / max(np.finfo(float).eps, denom)

    res, _ = integrate.quad(integrand, 0.0, R, limit=200)
    return (2.0 / (R**2)) * res


def compute_rho_eq14(E1, u, v, alpha=0.0, eps=1e-10):
    """Roze (2023) eq. (14): ρ ≈ 1 + (ε₁/(1 − u ε₁)) ((u + v + α)/2)."""
    denom = 1.0 - u * E1
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    return 1.0 + (E1 / denom) * ((u + v + alpha) / 2.0)


# ----- Fixed-width ASCII table (single-byte pipes; aligns in all terminals)
W_LAB = 20  # fits "Theory (neg. bin.)"
W_NUM = 14


def _hline():
    parts = ["-" * W_LAB] + ["-" * W_NUM] * 5
    return "  +" + "+".join(parts) + "+"


def _row_cells(cells):
    """cells: 6 strings, each already width-padded."""
    return "  | " + " | ".join(cells) + " |"


def _pad_num(x, decimals, empty="—"):
    if np.isfinite(x):
        return f"{x:>{W_NUM}.{decimals}f}"
    return f"{empty:>{W_NUM}}"


# ----- Per-row tables: mean, variance, ρ, skewness, excess kurtosis
n_rows = len(df)

for i in range(n_rows):
    Ri = float(R[i])
    e1 = compute_E1(Ri, U_ROZE)
    roze_mean = (U_ROZE - V_ROZE - ALPHA_ROZE) / BETA_ROZE if BETA_ROZE > 0 else np.nan
    # Var(n) ≈ n̅ + E₁ ((u+v+α)/2) n̅  from (5) and (13)
    roze_var = (
        roze_mean * (1.0 + e1 * (U_ROZE + V_ROZE + ALPHA_ROZE) / 2.0)
        if np.isfinite(roze_mean) and roze_mean > 0
        else np.nan
    )
    roze_rho = roze_var / roze_mean if np.isfinite(roze_mean) and roze_mean > 0 else np.nan
    rho_eq14 = compute_rho_eq14(e1, U_ROZE, V_ROZE, ALPHA_ROZE)

    print(f"\n[{i + 1}/{n_rows}]  R = {Ri:g}"
          f"   (approx. / neg. bin. rows at the {LIFE_STAGE}-transposition/excision stage)")
    print(_hline())
    print(
        _row_cells(
            [
                f"{'':{W_LAB}}",
                f"{'Mean':^{W_NUM}}",
                f"{'Var.':^{W_NUM}}",
                f"{'ρ':^{W_NUM}}",
                f"{'skew.':^{W_NUM}}",
                f"{'Ex. kurt.':^{W_NUM}}",
            ]
        )
    )
    print(_hline())

    sim_rho = sim_var_st[i] / sim_mean_st[i] if sim_mean_st[i] > 0 else np.nan
    app_rho = app_var_st[i] / app_mean_st[i] if app_mean_st[i] > 0 else np.nan
    nb_rho = nb_var_st[i] / nb_mean_st[i] if nb_mean_st[i] > 0 else np.nan

    if show_sim_shape:
        print(
            _row_cells(
                [
                    f"{'stochastic sim.':<{W_LAB}}",
                    _pad_num(sim_mean_st[i], 5),
                    _pad_num(sim_var_st[i], 5),
                    _pad_num(sim_rho, 5),
                    _pad_num(sim_skew_st[i], 4),
                    _pad_num(sim_exkurt_st[i], 4),
                ]
            )
        )
    else:
        print(
            _row_cells(
                [
                    f"{'stochastic sim.':<{W_LAB}}",
                    _pad_num(sim_mean_st[i], 5),
                    _pad_num(sim_var_st[i], 5),
                    _pad_num(sim_rho, 5),
                    f"{'—':>{W_NUM}}",
                    f"{'—':>{W_NUM}}",
                ]
            )
        )
    if show_nb_shape:
        print(
            _row_cells(
                [
                    f"{'Theory (neg. bin.)':<{W_LAB}}",
                    _pad_num(nb_mean_st[i], 5),
                    _pad_num(nb_var_st[i], 5),
                    _pad_num(nb_rho, 5),
                    _pad_num(nb_skew_st[i], 4),
                    _pad_num(nb_exkurt_st[i], 4),
                ]
            )
        )
    else:
        print(
            _row_cells(
                [
                    f"{'Theory (neg. bin.)':<{W_LAB}}",
                    _pad_num(nb_mean_st[i], 5),
                    _pad_num(nb_var_st[i], 5),
                    _pad_num(nb_rho, 5),
                    f"{'—':>{W_NUM}}",
                    f"{'—':>{W_NUM}}",
                ]
            )
        )
    print(
        _row_cells(
            [
                f"{'theory (approx.)':<{W_LAB}}",
                _pad_num(app_mean_st[i], 5),
                _pad_num(app_var_st[i], 5),
                _pad_num(app_rho, 5),
                f"{'—':>{W_NUM}}",
                f"{'—':>{W_NUM}}",
            ]
        )
    )
    print(
        _row_cells(
            [
                f"{'Roze (2023)':<{W_LAB}}",
                _pad_num(roze_mean, 5),
                _pad_num(roze_var, 5),
                _pad_num(roze_rho, 5),
                f"{'—':>{W_NUM}}",
                f"{'—':>{W_NUM}}",
            ]
        )
    )
    print(_hline())
    print(f"  ε₁ (E₁) at R = {Ri:g}  →  {e1:.6g}   |   ρ_eq14 = {rho_eq14:.6g}")

if not has_skew:
    print("\nNote: skewness / excess kurtosis columns not in CSV.")
