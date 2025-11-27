import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# CSV file paths - modify these to your actual file paths
csv_file_path_1 = "path/s:0.0001_N:10000_r:1_beta:0.032_delta:0.0064_exp.csv"
csv_file_path_2 = "path/s:0.0001_N:10000_r:1_beta:0.1_delta:0.07_exp.csv"


# Load both datasets
data1 = pd.read_csv(csv_file_path_1)
data2 = pd.read_csv(csv_file_path_2)

# Define distinct styles for each plot
# Plot 1 - Exponential fitness function
mean_styles_1 = {
    'Simulated Mean': {'color': '#0173B2', 'linestyle': '-', 'linewidth': 5.5},                    # Bright blue, solid
    'NB Mean': {'color': 'red', 'linestyle': '-', 'linewidth': 5.0},                      # Green, dash-dot
}

variance_styles_1 = variance_styles = {
    'Simulated Variance': {'color': '#4B0082', 'linestyle': '-', 'linewidth': 1.8},               # Web-safe red, solid
    'NB Variance': {'color': '#000000', 'linestyle': '-', 'linewidth': 6.5},                 # Deep pink, dash-dot
    # 'Charles Mean': {'color': '#FF4500', 'linestyle': ':', 'linewidth': 7.0}      # Orange red, dash-dot-dot pattern
    'Charles Mean': {'color': '#299438', 'linestyle': ':', 'linewidth': 7.0}
}

# Plot 2 styles - Synergistic fitness function
mean_styles_2 = {
    'Simulated Mean': {'color': '#0173B2', 'linestyle': '-', 'linewidth': 5.5},                    # Bright blue, solid
    'NB Mean': {'color': 'red', 'linestyle': '-', 'linewidth': 5.5},                      # Green, dash-dot
}

variance_styles_2 = {
    'Simulated Variance': {'color': '#4B0082', 'linestyle': '-', 'linewidth': 1.8},               # Web-safe red, solid
    'NB Variance': {'color': '#000000', 'linestyle': '-', 'linewidth': 6.5},                 # Deep pink, dash-dot
    # 'Charles Mean': {'color': '#FF4500', 'linestyle': ':', 'linewidth': 6.5}      # Orange red, dash-dot-dot pattern
    'Charles Mean': {'color': '#299438', 'linestyle': ':', 'linewidth': 6.5}
}


def plot_data_set1(data, ax, title, mean_styles, variance_styles):
    """Function to plot exponential fitness function"""
    # Extract columns
    generation = data['generation']
    mean_TE_per_individual = data['Mean_X_Count']
    NB_Mean = data['NB_Mean']
    Simulated_Variance = data['Sim_Variance']
    NB_Variance = data['NB_Variance']
    Charles_Mean = data['Charles_Mean']

    ####### To calculate the final values for the simulated mean and variance (exponential)
    t0 = 500
    k = 100

    indices1 = np.arange(t0, len(mean_TE_per_individual), k)
    indices2 = np.arange(t0, len(Simulated_Variance), k)

    mean_TE_piecewise = mean_TE_per_individual[indices1]
    sim_variance_piecewise = Simulated_Variance[indices2]

    final_average_mean_TE = mean_TE_piecewise.mean()
    final_average_simulated_variance = sim_variance_piecewise.mean()

    # Plot simulated vs approximate comparisons
    ax.plot(generation[:2000], NB_Mean[:2000], label='Theoretical Mean', **mean_styles['NB Mean'])
    ax.plot(generation[:2000], mean_TE_per_individual[:2000], label='Simulated Mean', **mean_styles['Simulated Mean'], alpha=0.9)
    ax.plot(generation[:2000], Simulated_Variance[:2000], label='Simulated Variance', **variance_styles['Simulated Variance'], alpha=0.9)
    ax.plot(generation[:2000], NB_Variance[:2000], label='Theoretical Variance', **variance_styles['NB Variance'])
    ax.plot(generation[:2000], Charles_Mean[:2000], label='Charlesworth (1983) Mean', **variance_styles['Charles Mean'])

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
    #           prop={'weight':'bold'})

    # Print final values for this dataset
    print(f"\n{title}:")
    print(f"Simulated Mean: {mean_TE_per_individual.iloc[-1]:.6f}")
    print(f"NB Mean: {NB_Mean.iloc[-1]:.6f}")
    print(f"Simulated Variance: {Simulated_Variance.iloc[-1]:.6f}")
    print(f"NB Variance: {NB_Variance.iloc[-1]:.6f}")
    print(f"Charles Mean: {Charles_Mean.iloc[-1]:.6f}")

    print(f"\nSimulated Mean_1 (piecewise_average): {final_average_mean_TE:.2f}")
    print(f"Simulated Variance_1 (piecewise_average): {final_average_simulated_variance:.2f}")


def plot_data_set2(data, ax, title, mean_styles, variance_styles):
    """Function to plot the synergistic fitness function"""
    # Extract columns
    generation = data['generation']
    mean_TE_per_individual = data['Mean_X_Count']
    NB_Mean = data['NB_Mean']
    Simulated_Variance = data['Sim_Variance']
    NB_Variance = data['NB_Variance']
    Charles_Mean = data['Charles_Mean']

    ####### To calculate the final values for the simulated mean and variance (exponential)
    t0 = 500
    k = 100

    indices1 = np.arange(t0, len(mean_TE_per_individual), k)
    indices2 = np.arange(t0, len(Simulated_Variance), k)

    mean_TE_piecewise = mean_TE_per_individual[indices1]
    sim_variance_piecewise = Simulated_Variance[indices2]

    final_average_mean_TE = mean_TE_piecewise.mean()
    final_average_simulated_variance = sim_variance_piecewise.mean()

    # Plot exact vs approximate comparisons
    ax.plot(generation[:2000], NB_Mean[:2000], label='Theoretical Mean', **mean_styles['NB Mean'])
    ax.plot(generation[:2000], mean_TE_per_individual[:2000], label='Simulated Mean', **mean_styles['Simulated Mean'], alpha=0.9)
    ax.plot(generation[:2000], Simulated_Variance[:2000], label='Simulated Variance', **variance_styles['Simulated Variance'],
            alpha=0.85)
    ax.plot(generation[:2000], NB_Variance[:2000], label='Theoretical Variance', **variance_styles['NB Variance'])
    ax.plot(generation[:2000], Charles_Mean[:2000], label='Charlesworth (1983) Mean', **variance_styles['Charles Mean'])

    # Formatting
    ax.set_xlabel("Generation", fontsize=22, fontweight='bold')
    # ax.set_ylabel("Mean Copies per Individual", fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=20, fontweight='bold', pad=5)
    ax.grid(True, linestyle='--', alpha=0, linewidth=0.5)

    # Make axis numbers bold
    ax.tick_params(axis='both', which='major', labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # # Add legend
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
    # print(f"Simulated Mean: {mean_TE_per_individual.iloc[-1]:.6f}")
    # print(f"NB Mean: {NB_Mean.iloc[-1]:.6f}")
    # print(f"NB Variance: {NB_Variance.iloc[-1]:.6f}")
    print(f"Simulated Mean: {mean_TE_per_individual.iloc[-1]:.6f}")
    print(f"NB Mean: {NB_Mean.iloc[-1]:.6f}")
    print(f"Simulated Variance: {Simulated_Variance.iloc[-1]:.6f}")
    print(f"NB Variance: {NB_Variance.iloc[-1]:.6f}")
    print(f"Charles Mean: {Charles_Mean.iloc[-1]:.6f}")

    print(f"\nSimulated Mean_2 (piecewise_average): {final_average_mean_TE:.2f}")
    print(f"Simulated Variance_2 (piecewise_average): {final_average_simulated_variance:.2f}")


# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7), dpi=120)

# Plot first dataset - Exponential fitness function
plot_data_set1(data1, ax1, r"$\bf{u = 0.032}$", mean_styles_1, variance_styles_1)

# Plot second dataset - Synergistic fitness function
plot_data_set2(data2, ax2, r"$\bf{u = 0.10}$", mean_styles_2, variance_styles_2)


# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()


