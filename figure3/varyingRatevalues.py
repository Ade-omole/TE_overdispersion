"""Print summary tables for the u-sweep (simulation vs theory)."""

from pathlib import Path

import numpy as np
import pandas as pd

TE_ROOT = Path(__file__).resolve().parent.parent
CSV_DIR = TE_ROOT / "csv_files"

R_val = 10.0
fitness = "exp"
CSV_BASENAME = "sweep_u_more_values_vConst:0.0001_6_binomial"

csv_path = CSV_DIR / f"{CSV_BASENAME}_{fitness}.csv"
if not csv_path.is_file():
    print(f"CSV not found: {csv_path}")
    print("Run varyingRatesim.jl first, or change CSV_BASENAME / fitness.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
print(f"Loaded {csv_path}, shape {df.shape}")

u_vals = df["u"].values
v_vals = df["v"].values

# ----- Post-transposition (as measured in the simulation / CSV)
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

# ----- Pre-transposition simulation (the "_m0" columns)
has_m0 = "sim_μ_m0" in df.columns
if has_m0:
    sim_mean_m0 = df["sim_μ_m0"].values
    sim_var_m0 = df["sim_σ²_m0"].values
    sim_p_m0 = df["sim_p_m0"].values

# ----- Theory (Neg. bin.), pre-transposition
S_VAL = 1e-4  # matches fit="exp", w(n) = exp(-s n^2)


def pretransp_nb_theory(mu, sig2, s=S_VAL):
    """Pre-transposition NB theory (Eqs. 6-7): returns (mean, variance, rho)."""
    p = mu / sig2
    sigma = np.sqrt(sig2)
    gamma = ((2.0 - p) / p) / sigma
    kappa = (6.0 * (1.0 - p) + p**2) / sig2 + 3.0

    beta1, beta2 = -2.0 * s * mu, -2.0 * s
    e_m0 = mu + beta1 * sig2 + 0.5 * beta2 * gamma * sig2**1.5
    var_m0 = (
        0.5 * (mu + sig2)
        + 0.5 * beta1 * (sig2 + gamma * sig2**1.5)
        + 0.25 * beta2 * (gamma * sig2**1.5 + (kappa - 1.0) * sig2**2)
    )
    return e_m0, var_m0, var_m0 / e_m0


# ----- Theory (approx.), pre-transposition
def approx_pretransp_rho(u, v):
    """rho^app_{m'} = (1 + 2v) / (1 - 2u + 2v)."""
    denom = 1.0 - 2.0 * u + 2.0 * v
    return (1.0 + 2.0 * v) / denom if denom != 0 else np.nan


# ----- Roze (2023): mean and variance from the paper's own equations
ALPHA_ROZE = 0.0
BETA_ROZE = 2.0 * S_VAL  # matches fit="exp" (w(n) = exp(-s n^2) <=> W = exp(-alpha n - beta n^2/2))
roze_mean = (u_vals - v_vals - ALPHA_ROZE) / BETA_ROZE

# Pre-transposition variance: Roze (2023) eq. (5), Var(n) ~= n̄ + 2*sum(D_ij),
roze_var_preT = roze_mean * (1.0 + 2.0 * u_vals) - BETA_ROZE * roze_mean**2

# Post-transposition variance: same mean (eq. 11) + review-file eq. (7),
roze_var_postT = roze_mean * (1.0 + 4.0 * u_vals) - BETA_ROZE * roze_mean**2

# rho columns: computed directly from the mean/variance above (rather than the
roze_rho_preT = roze_var_preT / roze_mean
roze_rho_postT = roze_var_postT / roze_mean

n_total = len(u_vals)
w1 = 28
header = f"  ┌{'─' * w1}┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐"
colhead = f"  │{'':^{w1}}│     Mean     │   Variance   │      ρ       │   Skewness   │  Ex. Kurtosis│"
midrule = f"  ├{'─' * w1}┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤"
footer = f"  └{'─' * w1}┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘"


def fmt_row(label, mean, var, rho, skew=None, exkurt=None):
    mean_s = f"{mean:12.2f}" if mean is not None and np.isfinite(mean) else f"{'—':^12}"
    var_s = f"{var:12.2f}" if var is not None and np.isfinite(var) else f"{'—':^12}"
    rho_s = f"{rho:12.4f}" if rho is not None and np.isfinite(rho) else f"{'—':^12}"
    skew_s = f"{skew:12.2f}" if skew is not None and np.isfinite(skew) else f"{'—':^12}"
    exkurt_s = f"{exkurt:12.2f}" if exkurt is not None and np.isfinite(exkurt) else f"{'—':^12}"
    return f"  │{label:<{w1}}│ {mean_s} │ {var_s} │ {rho_s} │ {skew_s} │ {exkurt_s} │"


for i in range(n_total):
    u, v = u_vals[i], v_vals[i]

    # Shared skewness / excess kurtosis -- same values used in both the
    sk = sim_skew[i] if has_skew else None
    ek = sim_exkurt[i] if has_skew else None
    nsk = nb_skew[i] if has_skew else None
    nek = nb_exkurt[i] if has_skew else None

    # ----- Pre-transposition table
    print(f"\n[{i + 1}/{n_total}]  u = {u}  (v = {v}, fitness = {fitness}) — Pre-transposition")
    print(header)
    print(colhead)
    print(midrule)
    if has_m0:
        sm0, sv0, sp0 = sim_mean_m0[i], sim_var_m0[i], sim_p_m0[i]
        sim_rho0 = sv0 / sm0 if sm0 > 0 else np.nan
        print(fmt_row("Simulation- Pre T.", sm0, sv0, sim_rho0, sk, ek))
    else:
        print(fmt_row("Simulation- Pre T.", None, None, None))

    e_m0, var_m0, rho_nb_preT = pretransp_nb_theory(nb_mean[i], nb_var[i])
    print(fmt_row("Theory (Neg. bin.)- Pre T.", e_m0, var_m0, rho_nb_preT, nsk, nek))

    # Theory (approx.) only gives a formula for rho (the dispersion ratio);
    rho_app_preT = approx_pretransp_rho(u, v)
    print(fmt_row("Theory (approx.)- Pre T.", app_mean[i], app_mean[i] * rho_app_preT, rho_app_preT))

    print(fmt_row("Roze (2023)- Pre T.", roze_mean[i], roze_var_preT[i], roze_rho_preT[i]))
    print(footer)

    # ----- Post-transposition table
    print(f"\n[{i + 1}/{n_total}]  u = {u}  (v = {v}, fitness = {fitness}) — Post-transposition")
    print(header)
    print(colhead)
    print(midrule)
    sm, sv, sp = sim_mean[i], sim_var[i], sim_p[i]
    sim_rho = sv / sm if sm > 0 else np.nan
    print(fmt_row("Simulation- Post T.", sm, sv, sim_rho, sk, ek))

    nm, nv = nb_mean[i], nb_var[i]
    nb_rho = nv / nm if nm > 0 else np.nan
    print(fmt_row("Theory (Neg. bin.)- Post T.", nm, nv, nb_rho, nsk, nek))

    am, av = app_mean[i], app_var[i]
    app_rho = av / am if am > 0 else np.nan
    print(fmt_row("Theory (approx.)- Post T.", am, av, app_rho))

    print(fmt_row("Roze (2023)- Post T.", roze_mean[i], roze_var_postT[i], roze_rho_postT[i]))
    print(footer)

if not has_skew:
    print("\nNote: skewness / excess kurtosis columns not in CSV. Re-run varyingRatesim.jl if needed.")
