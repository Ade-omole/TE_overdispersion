"""Print summary tables for the varying-u sweep (simulation vs theory)."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate

TE_ROOT = Path(__file__).resolve().parent.parent
CSV_DIR = TE_ROOT / "csv_files"

R_val = 10.0
fitness = "exp"
CSV_BASENAME = "sweep_u_roze_more_values_vConst:0.0001_5_binomial"

csv_path = CSV_DIR / f"{CSV_BASENAME}_{fitness}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run varyingRatesim.jl first, or change CSV_BASENAME / fitness.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

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

# Roze (2023): per-row u, v; table rho from eqs. (5)+(13); footer uses eq. (14).
ALPHA_ROZE = 0.0
BETA_ROZE = 2.0 / 10**4


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
    """Roze (2023) eq. (14): ρ ≈ 1 + (ε₁/(1 − u ε₁)) ((u + v + α)/2)."""
    denom = 1.0 - u * E1
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    return 1.0 + E1 * ((u + v + alpha) / 2.0)


u_vals = df["u"].values
v_vals = df["v"].values
E1_vals = np.array([compute_E1(R_val, u) for u in u_vals])
rho_roze = np.array([compute_rho_eq14(E1, u, v) for E1, u, v in zip(E1_vals, u_vals, v_vals)])

roze_mean_arr = np.array(
    [
        (u - v - ALPHA_ROZE) / BETA_ROZE if BETA_ROZE > 0 else np.nan
        for u, v in zip(u_vals, v_vals)
    ]
)
roze_var_arr = np.array(
    [
        m * (1.0 + e1 * (u + v + ALPHA_ROZE) / 2.0)
        if np.isfinite(m) and m > 0
        else np.nan
        for m, e1, u, v in zip(roze_mean_arr, E1_vals, u_vals, v_vals)
    ]
)
roze_rho_vm_arr = np.where(
    np.isfinite(roze_mean_arr) & (roze_mean_arr > 0),
    roze_var_arr / roze_mean_arr,
    np.nan,
)

n_total = len(u_vals)
for i in range(n_total):
    u, v = u_vals[i], v_vals[i]
    e1 = E1_vals[i]
    rm, rv, rvm = roze_mean_arr[i], roze_var_arr[i], roze_rho_vm_arr[i]
    r14 = rho_roze[i]
    inv_sim = 1.0 / max(sim_p[i], 1e-12)
    inv_app = 1.0 / max(app_p[i], 1e-12)
    inv_nb = 1.0 / max(nb_p[i], 1e-12)

    w1 = 22
    print(f"\n[{i + 1}/{n_total}]  u = {u}  (v = {v}, fitness = {fitness}, R = {R_val}) [Roze]")
    print("  ┌──────────────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐")
    print(f"  │{'':^{w1}}│     Mean     │   Variance   │      ρ       │   Skewness   │  Ex. Kurtosis│")
    print("  ├──────────────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤")
    sm, sv, sp = sim_mean[i], sim_var[i], sim_p[i]
    sim_rho = sv / sm if sm > 0 else np.nan
    if has_skew:
        print(
            f"  │{'Simulation':<{w1}}│ {sm:12.2f} │ {sv:12.2f} │ {sim_rho:12.4f} │ {sim_skew[i]:12.2f} │ {sim_exkurt[i]:12.2f} │"
        )
    else:
        print(f"  │{'Simulation':<{w1}}│ {sm:12.2f} │ {sv:12.2f} │ {sim_rho:12.4f} │      —       │      —       │")
    nm, nv, np_ = nb_mean[i], nb_var[i], nb_p[i]
    nb_rho = nv / nm if nm > 0 else np.nan
    if has_skew:
        print(
            f"  │{'Theory (Neg. bin.)':<{w1}}│ {nm:12.2f} │ {nv:12.2f} │ {nb_rho:12.4f} │ {nb_skew[i]:12.2f} │ {nb_exkurt[i]:12.2f} │"
        )
    else:
        print(f"  │{'Theory (Neg. bin.)':<{w1}}│ {nm:12.2f} │ {nv:12.2f} │ {nb_rho:12.4f} │      —       │      —       │")
    am, av, ap = app_mean[i], app_var[i], app_p[i]
    app_rho = av / am if am > 0 else np.nan
    print(f"  │{'Theory (approx.)':<{w1}}│ {am:12.2f} │ {av:12.2f} │ {app_rho:12.4f} │      —       │      —       │")
    print(
        f"  │{'Roze (2023)':<{w1}}│ {rm:12.2f} │ {rv:12.2f} │ {rvm:12.4f} │      —       │      —       │"
    )
    print("  └──────────────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘")
    print(
        f"  1/p (ρ⁻¹):  sim = {inv_sim:.6g}   Theory (approx.) = {inv_app:.6g}   "
        f"Theory (Neg. bin.) = {inv_nb:.6g}   Roze (2023, eq.5+13) = {rvm:.6g}"
    )
    print(
        f"  ε₁ (E₁) at R={R_val}, u={u:g} → {e1:.6g}   |   rho_combine_eq14 (Roze eq. 14) = {r14:.6g}   "
        f"(differs from Roze table ρ = var/mean when eq. 14 ≠ eq. 5+13 closure)"
    )

if not has_skew:
    print("\nNote: skewness / excess kurtosis columns not in CSV. Re-run varyingRatesim.jl if needed.")
