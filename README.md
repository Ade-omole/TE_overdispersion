# Population genetics of active transposable elements: overdispersion arises naturally from transposition-deletion-selection balance

Code and data to reproduce the stochastic simulations, theoretical closures, and empirical analyses. Each subfolder corresponds to one main-text or appendix figure; reference PDFs are included alongside the scripts.

## Repository layout

```
TE_overdispersion/
├── csv_files/          # Simulation summaries (shared across figure folders)
├── figure1/            # Fig. 1 — approximate & NB closure time series (appNB.pdf)
├── figure2/            # Fig. 2 — skewness & excess kurtosis (skew_and_exkurt.pdf)
├── figure3/            # Fig. 3 — transposition-rate sweep (varyingRateFig.pdf)
├── figure4/            # Fig. 4 — selection-strength sweep (varyingSelection.pdf)
├── figure5/            # Figs. 5 & 6 — DPGP3 data (mean–variance–ρ and copy-number panels)
├── figureD1/           # Appendix Fig. D.1 — map length R sweep (varyingRFig.pdf)
└── figureE1/           # Appendix Fig. E.1 — fitness curvature sweep (varyingFitness_exp.pdf)
```

Scripts read/write paths relative to this folder: Julia simulations write CSVs to `csv_files/`; plotting scripts save PDFs into their own figure folder.

## Life-cycle stage

TE copy number is measured at two points in the life cycle, and the two differ:

- **post-transposition/excision** — after transposition and excision have acted, the standing population distribution (main text `n`, dispersion ratio `ρ ≈ 1 + 3u`);
- **pre-transposition/excision** — immediately after recombination, before the current generation's transposition and excision (main text `m'`, `ρ_m' ≈ 1 + u`), the stage to which the Roze (2023) prediction applies.

The Julia sweeps record both in a single run: post-transposition columns (`sim_μ`, `sim_σ²`, `sim_p`, …) and pre-transposition columns suffixed `_m0` (`sim_μ_m0`, `sim_σ²_m0`, `sim_p_m0`, `sim_skew_m0`, `sim_exkurt_m0`). `figureD1/varyingRFig.py` and `figureD1/varyingRvalues.py` expose a `LIFE_STAGE` switch (`"post"` or `"pre"`) that selects the stage for both theory curves and simulation markers.

Most plotting scripts also carry a `NB_THEORY_SOURCE` (or `THEORY_SOURCE`) switch choosing between the numerically integrated closure and the closed-form equilibria of the main text — Eq. (14) for `ρ̂_nb`, Eq. (15) for mean and variance, Eq. (16) for skewness and excess kurtosis, and Eqs. (23)–(24) for the pre-transposition moments.

## Requirements

**Julia** (1.6+ recommended) with:

```julia
using Pkg
Pkg.add(["CSV", "DataFrames", "Distributions", "QuadGK", "Statistics"])
```

(`Random`, `Printf`, and `Statistics` are in the standard library.)

**Python** (3.8+ recommended):

```bash
pip install numpy pandas matplotlib scipy jupyter
```

## Quick start (figures only)

Pre-computed CSVs are included in `csv_files/`, so you can regenerate figures without re-running long simulations:

```bash
cd TE_overdispersion

python figure1/appNBfig.py        # Fig. 1
python figure2/skew_exKurtFig.py  # Fig. 2 (same CSV as Fig. 1)
python figure3/varyingRateFig.py  # Fig. 3
python figure4/varyingSelection.py # Fig. 4
python figureD1/varyingRFig.py    # Fig. D.1
python figureE1/varyingFitness.py # Fig. E.1
```

For **Figs. 5–6**, open and run all cells in `figure5/DPGP3_TE_matrix_preview.ipynb` (from `figure5/` or the repo root).

Each command overwrites the corresponding PDF in the figure folder. Compare with the bundled reference PDFs to check reproducibility.

## Full reproduction (simulations + figures)

Run Julia simulations first where noted, then the Python plotting scripts. Simulations can take hours depending on hardware.

| Figure | Step 1 — simulation | Step 2 — plot / analysis | Output PDF |
|--------|---------------------|--------------------------|------------|
| 1 | `julia figure1/appNB.jl` | `python figure1/appNBfig.py` | `figure1/appNB.pdf` |
| 2 | *(uses Fig. 1 CSV)* | `python figure2/skew_exKurtFig.py` | `figure2/skew_and_exkurt.pdf` |
| 3 | `julia figure3/varyingRatesim.jl` | `python figure3/varyingRateFig.py` | `figure3/varyingRateFig.pdf` |
| 4 | `julia figure4/varyingSelection.jl` | `python figure4/varyingSelection.py` | `figure4/varyingSelection.pdf`, `figure4/varyingSelection_table.pdf` |
| 5–6 | — | Run `figure5/DPGP3_TE_matrix_preview.ipynb` | `figure5/DPGP3_TE_mean_variance_rho.pdf`, `figure5/DPGP3_TE_copynumber_POGO_JOCKEY_MDG1.pdf` |
| D.1 | `julia figureD1/varyingRsim.jl` **and** `julia figureD1/sweep_R_ectFit.jl` | `python figureD1/varyingRFig.py` | `figureD1/varyingRFig.pdf` |
| E.1 | `julia figureE1/varyingFitnesssim.jl` | `python figureE1/varyingFitness.py` | `figureE1/varyingFitness_exp.pdf` |

Fig. D.1 needs both Julia sweeps: `varyingRsim.jl` for the Gaussian fitness `w(n) = exp(-s n²)` and `sweep_R_ectFit.jl` for the Roze (2023) ectopic-pair fitness `w = exp(-β n_p)`. `varyingFitnesssim.jl` writes one CSV per fitness-curvature value `γ`.

### Optional: equilibrium tables (supplementary)

These print numeric summaries to the terminal (no PDF):

```bash
python figure3/varyingRatevalues.py   # u-sweep equilibrium moments
python figureD1/varyingRvalues.py     # R-sweep equilibrium moments
```

## Bundled data files

| File | Used by |
|------|---------|
| `csv_files/appNB_s:…_R:free_exp.csv` | Figs. 1 & 2 |
| `csv_files/sweep_u_more_values_vConst:0.0001_6_binomial_exp.csv` | Fig. 3 (`varyingRateFig.py`, `varyingRatevalues.py`) |
| `csv_files/varyingSelection_binomial_exp_u0.01_v0.0001_4.csv` | Fig. 4 |
| `csv_files/sweep_R_topright_burnin_binomial_pre_exp.csv` | Fig. D.1 (Gaussian fitness simulations) |
| `csv_files/sweep_R_ectFit_burnin_binomial_pre_exp.csv` | Fig. D.1 (Roze ectopic-pair fitness simulations) |
| `csv_files/varyingFitness_sim_exp_gamma{2.5,4.0,7.0,10.0}.csv` | Fig. E.1 |
| `figure5/DPGP3_TE_insertion_depletion_summary_everything_masked_noHet_bool.txt` | Figs. 5 & 6 (Lee et al. 2022 DPGP3 data) |

## Model summary

- **Stochastic simulations:** diploid Wright–Fisher–Moran-style updates with binomial transposition/excision, Gaussian fitness `w(n) = exp(-s n²)`, and free (binomial) recombination unless stated otherwise.
- **Theory:** approximate moment closure and negative-binomial parametric closure, available either as numerically integrated ODEs or as closed-form equilibria (see the switches above).
- **Empirical data:** per-strain TE insertion burdens from DPGP3 Zambian *D. melanogaster* (euchromatic insertions; see notebook for filtering).

## Citation

If you use this code or data, please cite:

Omole, A. D., & Czuppon, P. Population genetics of active transposable elements: overdispersion arises naturally from transposition-deletion-selection balance.

Empirical data source: [Lee et al. (2022)](https://academic.oup.com/genetics/article/220/2/iyab211/6458331), DPGP3 TE insertion matrix ([raw data](https://github.com/jumpingTE-LeeLab/TE_insertion_DPGP3_raw)).

Theory comparison: [Roze (2023)](https://academic.oup.com/genetics/article/224/2/iyad058/7109257?guestAccessKey=).
