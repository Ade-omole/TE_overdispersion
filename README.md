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
└── figureD1/           # Appendix Fig. D.1 — map length R sweep (varyingRFig.pdf)
```

Scripts read/write paths relative to this folder: Julia and Python simulations write CSVs to `csv_files/`; plotting scripts save PDFs into their own figure folder.

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

# Fig. 1
python figure1/appNBfig.py

# Fig. 2 (same CSV as Fig. 1)
python figure2/skew_exKurtFig.py

# Fig. 3
python figure3/varyingRateFig.py

# Fig. 4
python figure4/varyingSelection.py

# Fig. D.1
python figureD1/varyingRFig.py
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
| 4 | `julia figure4/varyingSelection.jl` | `python figure4/varyingSelection.py` | `figure4/varyingSelection.pdf` |
| 5–6 | — | Run `figure5/DPGP3_TE_matrix_preview.ipynb` | `figure5/DPGP3_TE_mean_variance_rho.pdf`, `figure5/DPGP3_TE_copynumber_POGO_JOCKEY_MDG1.pdf` |
| D.1 | `julia figureD1/varyingRsim.jl` | `python figureD1/varyingRFig.py` | `figureD1/varyingRFig.pdf` |

### Optional: equilibrium tables (supplementary)

These print numeric summaries to the terminal (no PDF):

```bash
python figure3/varyingRatevalues.py   # u-sweep equilibrium moments
python figureD1/varyingRvalues.py     # R-sweep equilibrium moments
```

`varyingSelection.py` can also write a supplementary-style table PDF (`varyingSelection_table.pdf`) when run with table output enabled in the script.

## Bundled data files

| File | Used by |
|------|---------|
| `csv_files/appNB_s:…_R:free_exp.csv` | Figs. 1 & 2 |
| `csv_files/sweep_u_roze_more_values_vConst:0.0001_4_binomial_exp.csv` | Fig. 3 (`varyingRateFig.py`) |
| `csv_files/varyingSelection_binomial_exp_u0.01_v0.0001_4.csv` | Fig. 4 |
| `csv_files/sweep_R_roze_topright_burnin_binomial_exp.csv` | Fig. D.1 (Gaussian fitness simulations) |
| `csv_files/sweep_R_ectFitRoze_burnin_binomial_exp.csv` | Fig. D.1 (Roze ectopic-pair fitness simulations) |
| `figure5/DPGP3_TE_insertion_depletion_summary_everything_masked_noHet_bool.txt` | Figs. 5 & 6 (Lee et al. 2022 DPGP3 data) |

## Model summary

- **Stochastic simulations:** diploid Wright–Fisher–Moran-style updates with binomial transposition/excision, Gaussian fitness `w(n) = exp(-s n²)`, and free (binomial) recombination unless stated otherwise.
- **Theory:** approximate moment closure and negative-binomial parametric closure (ODEs integrated alongside simulations).
- **Empirical data:** per-strain TE insertion burdens from DPGP3 Zambian *D. melanogaster* (euchromatic insertions; see notebook for filtering).

## Citation

If you use this code or data, please cite:

Omole, A. D., & Czuppon, P. A population genetics model explaining overdispersion in active transposable elements.

Empirical data source: [Lee et al. (2022)]([https://doi.org/10.1093/gbe/evac016](https://academic.oup.com/genetics/article/220/2/iyab211/6458331)), DPGP3 TE insertion matrix ([raw data](https://github.com/jumpingTE-LeeLab/TE_insertion_DPGP3_raw)).

Theory comparison: [Roze (2023)]([https://doi.org/10.1093/gbe/evad123](https://academic.oup.com/genetics/article/224/2/iyad058/7109257?guestAccessKey=)).
