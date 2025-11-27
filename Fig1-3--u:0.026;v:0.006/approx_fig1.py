import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# CSV file paths - modify these to your actual file paths
csv_file_path_1 = "path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_exp___.csv"
csv_file_path_2 = "path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_syn___.csv"

# Load both datasets
data1 = pd.read_csv(csv_file_path_1)
data2 = pd.read_csv(csv_file_path_2)

# Define distinct styles for each plot
# Plot 1 - Exponential fitness function
mean_styles_1 = {
    'Simulated Mean': {'color': '#0173B2', 'linestyle': '-', 'linewidth': 5.5},
    'Approx Mean': {'color': 'red', 'linestyle': '-', 'linewidth': 5.2},
}

variance_styles_1 = {
    'Simulated Variance': {'color': '#4B0082', 'linestyle': '-', 'linewidth': 1.9},
    'Approx Variance': {'color': '#000000', 'linestyle': '-', 'linewidth': 6.5},
    # 'Charles Mean': {'color': '#FF4500', 'linestyle': '--', 'linewidth': 6.5}
    'Charles Mean': {'color': '#299438', 'linestyle': '--', 'linewidth': 6.5}

}


# Plot 2 styles - Synergistic fitness function
mean_styles_2 = {
    'Simulated Mean': {'color': '#0173B2', 'linestyle': '-', 'linewidth': 5.5},
    'Approx Mean': {'color': 'red', 'linestyle': '-', 'linewidth': 5.3},
}

variance_styles_2 = {
    'Approx Variance': {'color': '#000000', 'linestyle': '-', 'linewidth': 6.5},
    # 'Charles Mean': {'color': '#FF4500', 'linestyle': ':', 'linewidth': 5.5},
    'Charles Mean': {'color': '#299438', 'linestyle': '--', 'linewidth': 7.0},
    'Simulated Variance': {'color': '#4B0082', 'linestyle': '-', 'linewidth': 1.9},
}


def plot_data_set1(data, ax, title, mean_styles, variance_styles):
    """Function to plot exponential fitness function"""
    # Extract columns
    generation = data['generation']
    mean_TE_per_individual = data['Mean_X_Count']
    Approx_Mean = data['Approx_Mean']
    Simulated_Variance = data['Sim_Variance']
    Approx_Variance = data['Approx_Variance']
    Charles_Mean = data['Charles_Mean']

    # Plot simulated vs approximate comparisons
    ax.plot(generation[:1000], Approx_Mean[:1000], label='Theoretical Mean', **mean_styles['Approx Mean'])
    ax.plot(generation[:1000], mean_TE_per_individual[:1000], label='Simulated Mean', **mean_styles['Simulated Mean'], alpha=0.9)
    ax.plot(generation[:1000], Simulated_Variance[:1000], label='Simulated Variance', **variance_styles['Simulated Variance'], alpha=0.85)
    ax.plot(generation[:1000], Approx_Variance[:1000], label='Theoretical Variance', **variance_styles['Approx Variance'])
    ax.plot(generation[:1000], Charles_Mean[:1000], label='Charlesworth (1983) Mean', **variance_styles['Charles Mean'])


    # Formatting
    ax.set_xlabel("Generation", fontsize=22, fontweight='bold')
    ax.set_ylabel("Copy number", fontsize=25, fontweight='bold')
    ax.set_title(title, fontsize=20, fontweight='bold', pad=5)
    ax.grid(True, linestyle='--', alpha=0, linewidth=0.5)

    # Make axis numbers bold
    ax.tick_params(axis='both', which='major', labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Add legend
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, labels,
    #           fontsize=10,
    #           frameon=True,
    #           fancybox=True,
    #           shadow=True,
    #           ncol=1,
    #           loc='best',
    #           prop={'size': 15, 'weight':'bold'})

    # Print final values for this dataset
    print(f"\n{title}:")
    print(f"Simulated Mean: {mean_TE_per_individual.iloc[-1]:.2f}")
    print(f"Approx Mean: {Approx_Mean.iloc[-1]:.2f}")
    print(f"Simulated Variance: {Simulated_Variance.iloc[-1]:.2f}")
    print(f"Approx Variance: {Approx_Variance.iloc[-1]:.2f}")
    print(f"Charles Mean: {Charles_Mean.iloc[-1]:.2f}")


def plot_data_set2(data, ax, title, mean_styles, variance_styles):
    """Function to plot the synergistic fitness function"""
    # Extract columns
    generation = data['generation']
    mean_TE_per_individual = data['Mean_X_Count']
    Approx_Mean = data['Approx_Mean']     ##### Mean from inverse dispersion p
    Simulated_Variance = data['Sim_Variance']
    Approx_Variance = data['Approx_Variance']   ##### Variance from inverse dispersion p
    Charles_Mean = data['Charles_Mean']

    # Plot exact vs approximate comparisons
    ax.plot(generation[:1000], Approx_Mean[:1000], label='Theoretical Mean', **mean_styles['Approx Mean'])
    ax.plot(generation[:1000], mean_TE_per_individual[:1000], label='Simulated Mean', **mean_styles['Simulated Mean'], alpha=0.9)
    ax.plot(generation[:1000], Approx_Variance[:1000], label='Theoretical Variance', **variance_styles['Approx Variance'])
    ax.plot(generation[:1000], Simulated_Variance[:1000], label='Simulated Variance', **variance_styles['Simulated Variance'], alpha=0.85)
    ax.plot(generation[:1000], Charles_Mean[:1000], label='Charlesworth (1983) Mean', **variance_styles['Charles Mean'])


    # Formatting
    ax.set_xlabel("Generation", fontsize=22, fontweight='bold')
    # ax.set_ylabel("Mean Copies per Individual", fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=20, fontweight='bold', pad=5)
    ax.grid(True, linestyle='--', alpha=0, linewidth=0.5)

    # Make axis numbers bold
    ax.tick_params(axis='both', which='major', labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Add legend
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, labels,
    #           fontsize=10,
    #           frameon=True,
    #           fancybox=True,
    #           shadow=True,
    #           ncol=1,
    #           loc='best',
    #           prop={'size': 14, 'weight':'bold'})

    # Print final values for this dataset
    print(f"\n{title}:")
    print(f"Simulated Mean: {mean_TE_per_individual.iloc[-1]:.2f}")
    print(f"Approx Mean: {Approx_Mean.iloc[-1]:.2f}")
    print(f"Approx Variance: {Approx_Variance.iloc[-1]:.2f}")

### For exponential
u = 0.026
v = 0.006
s = 1/10000
mean_app_theory = (u - v) * (1 - 2 * (u - v)) / (2 * s * (2 * u + 2 * v + 1))
var_app_theory = (u - v) / (2 * s)
print(f"Analytic Approx mean: {mean_app_theory:.2f}")
print(f"Analytic Approx variance: {var_app_theory:.2f}")
print(f"Ch. and Ch. mean and variance: {var_app_theory:.2f}")


# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7), dpi=120)

# Plot first dataset - Exponential fitness function
plot_data_set1(data1, ax1, r"$\bf{w(n) = exp(-sn^2)}$", mean_styles_1, variance_styles_1)

# Plot second dataset - Synergistic fitness function
plot_data_set2(data2, ax2, r"$\bf{w(n) = 1-sn^2}$", mean_styles_2, variance_styles_2)


# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()