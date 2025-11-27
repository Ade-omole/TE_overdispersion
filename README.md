# Transposable Element Overdispersion: Theory and Data Analysis

This repository contains the code and data for analyzing overdispersion in transposable element (TE) copy number distributions using Moran model simulations and empirical data analysis.

## Overview

This project investigates the statistical properties of transposable element copy number distributions in populations, focusing on overdispersion (variance exceeding the mean). The analysis combines:

- **Theoretical modeling**: Moran model simulations with selection and recombination
- **Approximate dynamics**: Comparison of approximate and negative binomial closure methods
- **Empirical data analysis**: Analysis of TE copy number data from *Drosophila melanogaster*

## Repository Structure

```
TE_overdispersion/
├── Fig1-3--u:0.026;v:0.006/          # Figures 1-3: Analysis for u=0.026, v=0.006
│   ├── simulation_u:0.026_exp.jl    # Julia simulation (exponential fitness)
│   ├── simulation_u:0.026_syn.jl    # Julia simulation (synergistic fitness)
│   ├── approx_fig1.py                # Figure 1: Approximate closure comparison
│   ├── nb_fig2.py                    # Figure 2: Negative binomial closure comparison
│   ├── skew&exKurt_fig3.py          # Figure 3: Skewness and excess kurtosis
│   ├── Fig1.pdf, Fig2.pdf, Fig3.pdf  # Generated figures
│   └── *.csv                         # Simulation output data
│
├── Fig4--u:0.032;u:0.1/              # Figure 4: Comparison of different u values
│   ├── simulation_u:0.032&0.1_exp.jl # Julia simulation
│   ├── fig4.py                       # Figure 4 visualization
│   ├── fig4.pdf                      # Generated figure
│   └── *.csv                         # Simulation output data
│
└── Fig5-6--Data/                     # Figures 5-6: Empirical data analysis
    ├── full_meta_data_mcgurk.py      # Metadata and statistics extraction
    ├── full_distribution_mcgurk_data.py # Distribution visualizations
    ├── sidebyside_Mcgurk_data_visualization.py # Side-by-side comparisons
    ├── mcgurk_data.xlsx              # McGurk et al. TE copy number data
    ├── Fig5.pdf, Fig6&supp.pdf       # Generated figures
```

## Requirements

### Julia Dependencies
- Julia 1.6 or higher
- Required packages (automatically installed):
  - `Random`
  - `Distributions`
  - `DataFrames`
  - `CSV`
  - `Statistics`

### Python Dependencies
- Python 3.7 or higher
- Required packages:
  - `pandas`
  - `matplotlib`
  - `numpy`
  - `scipy`
  - `seaborn` (for Fig5-6)
  - `openpyxl` (for reading Excel files)

Install Python dependencies:
```bash
pip install pandas matplotlib numpy scipy seaborn openpyxl
```

## Usage

### Running Simulations

#### Figure 1-3 Simulations

1. **Exponential fitness function**:
   ```bash
   julia Fig1-3--u:0.026;v:0.006/simulation_u:0.026_exp.jl
   ```
   This generates: `path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_exp___.csv`

2. **Synergistic fitness function**:
   ```bash
   julia Fig1-3--u:0.026;v:0.006/simulation_u:0.026_syn.jl
   ```
   This generates: `path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_syn___.csv`

**Note**: Update the `path/` placeholder in the Julia files to your desired output directory before running.

#### Figure 4 Simulation

```bash
julia Fig4--u:0.032;u:0.1/simulation_u:0.032&0.1_exp.jl
```

**Note**: This simulation uses u=0.032, v=0.0064. To generate data for u=0.1, v=0.07, modify the parameters in the file and run again.

### Generating Figures

#### Figures 1-3

1. **Figure 1 - Approximate closure comparison**:
   ```bash
   python Fig1-3--u:0.026;v:0.006/approx_fig1.py
   ```
   - Updates the `path/` placeholder to point to your CSV files
   - Compares approximate theoretical predictions with simulations for exponential and synergistic fitness functions

2. **Figure 2 - Negative binomial closure comparison**:
   ```bash
   python Fig1-3--u:0.026;v:0.006/nb_fig2.py
   ```
   - Updates the `path/` placeholder to point to your CSV files
   - Compares negative binomial theoretical predictions with simulations

3. **Figure 3 - Skewness and excess kurtosis**:
   ```bash
   python Fig1-3--u:0.026;v:0.006/skew&exKurt_fig3.py
   ```
   - Updates the `path/` placeholder to point to your CSV files
   - Visualizes skewness and excess kurtosis from simulations and theory

#### Figure 4

```bash
python Fig4--u:0.032;u:0.1/fig4.py
```
- Updates the `path/` placeholder to point to your CSV files
- Compares dynamics for different transposition rates (u=0.032 and u=0.1)

#### Figures 5-6

1. **Extract metadata and statistics**:
   ```bash
   python Fig5-6--Data/full_meta_data_mcgurk.py
   ```
   - Updates the `path/` placeholder to point to `mcgurk_data.xlsx`
   - Computes mean, variance, skewness, and excess kurtosis for each TE family

2. **Distribution visualizations**:
   ```bash
   python Fig5-6--Data/full_distribution_mcgurk_data.py
   ```
   - Updates the `path/` placeholder to point to `mcgurk_data.xlsx`
   - Generates histograms for all TE families

3. **Side-by-side comparisons**:
   ```bash
   python Fig5-6--Data/sidebyside_Mcgurk_data_visualization.py
   ```
   - Updates the `path/` placeholder to point to `mcgurk_data.xlsx`
   - Creates variance vs. mean scatter plots and dispersion index rankings

## Simulation Parameters

### Default Parameters (Fig1-3)
- **Population size (N)**: 10,000
- **Selection coefficient (s)**: 1/N = 0.0001
- **Recombination rate (r)**: 1 (free recombination)
- **Transposition rate (u)**: 0.026
- **Excision rate (v)**: 0.006
- **Chromosome length (L)**: 50,000 loci (exp) or 1,000 loci (syn)
- **Generations**: 20,000
- **Initial TE copies**: 5 per chromosome

### Figure 4 Parameters
- **Case 1**: u = 0.032, v = 0.0064
- **Case 2**: u = 0.10, v = 0.07
- Other parameters same as Fig1-3

## Data

### Simulation Output Format

CSV files contain the following columns:
- `generation`: Generation number
- `Mean_X_Count`: Simulated mean TE copy number
- `Approx_Mean`, `NB_Mean`: Theoretical mean predictions
- `Sim_Variance`: Simulated variance
- `Approx_Variance`, `NB_Variance`: Theoretical variance predictions
- `Charles_Mean`: Charlesworth (1983) mean prediction
- `Sim_Skewness`, `NB_Skewness`: Skewness values
- `Sim_Excess_kurtosis`, `NB_Excess_Kurtosis`: Excess kurtosis values
- Additional columns for intermediate calculations

### Empirical Data

The `mcgurk_data.xlsx` file contains TE copy number data from McGurk et al. (2021) for *Drosophila melanogaster*, including:
- DNA transposons: POGO, P-element, M4DM, BARI_DM, HOBO
- LTR transposons: ROO_LTR, DM412_LTR, NOMAD_LTR, MDG1_LTR, TRANSPAC_LTR, BURDOCK_LTR, TIRANT_LTR, TABOR_LTR
- Non-LTR transposons: Jockey, FW, I_DM, DOC, DOC6_DM

## Key Results

- **Overdispersion**: Variance exceeds mean for most TE families, indicating overdispersion
- **Fitness function effects**: Exponential and synergistic fitness functions show different dynamics
- **Theoretical predictions**: Negative binomial closure provides better agreement with simulations than approximate closure
- **Higher moments**: Skewness and excess kurtosis are well-predicted by negative binomial theory

## Citation

If you use this code or data, please cite:

Omole, A. D., & Czuppon, P. A population genetics model explaining overdispersion in active transposable elements.

