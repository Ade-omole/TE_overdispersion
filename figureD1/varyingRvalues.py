"""Print summary tables for the varying-R sweep (simulation vs theory)."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate

TE_ROOT = Path(__file__).resolve().parent.parent
CSV_DIR = TE_ROOT / "csv_files"

CSV_BASENAME = "sweep_R_roze_topright_burnin_binomial"
FITNESS = "exp"  # "exp" or "syn"

csv_path = CSV_DIR / f"{CSV_BASENAME}_{FITNESS}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run varyingRsim.jl first, or change CSV_BASENAME / FITNESS.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

R = df["R"].values

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

# Roze (2023): analytical rows use eqs. (11), (5)+(13), and (14).
ALPHA_ROZE = 0.0
U_ROZE = 0.01
V_ROZE = U_ROZE / 10.0
BETA_ROZE = U_ROZE / 10.0


def compute_E1(R, u, eps=1e-12):
    """E₁ (ε₁) for linear genetic map of length R Morgans (Roze appendix). R, u scalars."""
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
    """
    Roze (2023) eq. (14): rho ≈ 1 + (epsilon_1/(1 - u epsilon_1)) ((u + v + alpha)/2).
    """
    denom = 1.0 - u * E1
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    return 1.0 + (E1 / denom) * ((u + v + alpha) / 2.0)


W_LAB = 20
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

    print(f"\n[{i + 1}/{n_rows}]  R = {Ri:g}")
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

    sim_rho = sim_var[i] / sim_mean[i] if sim_mean[i] > 0 else np.nan
    app_rho = app_var[i] / app_mean[i] if app_mean[i] > 0 else np.nan
    nb_rho = nb_var[i] / nb_mean[i] if nb_mean[i] > 0 else np.nan

    if has_skew:
        print(
            _row_cells(
                [
                    f"{'stochastic sim.':<{W_LAB}}",
                    _pad_num(sim_mean[i], 5),
                    _pad_num(sim_var[i], 5),
                    _pad_num(sim_rho, 5),
                    _pad_num(sim_skew[i], 4),
                    _pad_num(sim_exkurt[i], 4),
                ]
            )
        )
    else:
        print(
            _row_cells(
                [
                    f"{'stochastic sim.':<{W_LAB}}",
                    _pad_num(sim_mean[i], 5),
                    _pad_num(sim_var[i], 5),
                    _pad_num(sim_rho, 5),
                    f"{'—':>{W_NUM}}",
                    f"{'—':>{W_NUM}}",
                ]
            )
        )
    if has_skew:
        print(
            _row_cells(
                [
                    f"{'Theory (neg. bin.)':<{W_LAB}}",
                    _pad_num(nb_mean[i], 5),
                    _pad_num(nb_var[i], 5),
                    _pad_num(nb_rho, 5),
                    _pad_num(nb_skew[i], 4),
                    _pad_num(nb_exkurt[i], 4),
                ]
            )
        )
    else:
        print(
            _row_cells(
                [
                    f"{'Theory (neg. bin.)':<{W_LAB}}",
                    _pad_num(nb_mean[i], 5),
                    _pad_num(nb_var[i], 5),
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
                _pad_num(app_mean[i], 5),
                _pad_num(app_var[i], 5),
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
