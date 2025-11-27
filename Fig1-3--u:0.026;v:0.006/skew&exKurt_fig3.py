import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# CSV file path
csv_file_path = "path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_syn___.csv"
# csv_file_path = "path/s:0.0001_N:10000_r:1_beta:0.026_delta:0.006_exp___.csv"

# Load the CSV file
try:
    df = pd.read_csv(csv_file_path)
    print("Data loaded successfully!")
    print(f"Shape: {df.shape}")
except FileNotFoundError:
    print(f"File not found: {csv_file_path}")
    print("Please check the file path and try again.")
    exit()

# Set up the plot style
plt.style.use('default')
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 12

# Remove first 40 data points and keep only 1000 time steps
df_plot = df.iloc[40:1040].reset_index(drop=True)

# Create x-axis (generation)
x = np.arange(len(df_plot))

# Define colors and line styles for better visibility
colors = ['#1f77b4', '#ff7f0e', '#000000']  # Blue, Orange, Black
linestyles = ['-', '-', '-']  # Solid, Solid, Solid
linewidths = [2.5, 6.5, 5.5]  # Different line widths

# Column name mapping for nicer legend + console output
label_map = {
    "Sim_Skewness": "Simulated Skewness",
    "NB_Skewness": "Theoretical (NB) Skewness",
    "Sim_Excess_kurtosis": "Simulated Excess Kurtosis",
    "NB_Excess_Kurtosis": "Theoretical (NB) Excess Kurtosis"
}

def plot_data_group(data, ax, columns, title, ylabel):
    """Function to plot a group of data on a given axis"""
    # Plot each column in the group
    for i, col in enumerate(columns):
        if col in data.columns:
            ax.plot(x, data[col],
                    color=colors[i % len(colors)],
                    linestyle=linestyles[i % len(linestyles)],
                    linewidth=linewidths[i % len(linewidths)],
                    alpha=0.9,
                    label=label_map.get(col, col),  # Use mapped label
                    marker='o' if len(data) < 50 else None,
                    markersize=4 if len(data) < 50 else 0)
        else:
            print(f"Warning: Column '{col}' not found in the dataset")

    # Formatting
    ax.set_xlabel('Generation', fontsize=22, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=22, fontweight='bold')
    ax.set_title(title, fontsize=20, fontweight='bold', pad=5)
    ax.grid(True, alpha=0)

    # Add legend inside the plot
    ax.legend(loc='best', prop={'size': 14, 'weight':'bold'})

    # Make axis numbers bold
    ax.tick_params(axis='both', which='major', labelsize=17, width=2)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Print final values for this group
    print(f"\n{title}:")
    for col in columns:
        if col in data.columns:
            final_value = data[col].iloc[-1]
            print(f"  {label_map.get(col, col)}: {final_value}")

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Define column groups
skewness_columns = ['Sim_Skewness', 'NB_Skewness']
kurtosis_columns = ['Sim_Excess_kurtosis', 'NB_Excess_Kurtosis']

# Plot skewness (left plot)
plot_data_group(df_plot, ax1, skewness_columns,
                "Skewness (NB)",
                "Skewness Value")

# Plot excess kurtosis (right plot)
plot_data_group(df_plot, ax2, kurtosis_columns,
                "Excess Kurtosis (NB)",
                "Excess Kurtosis Value")

####### To calculate the final values for the simulated skewness and excess kurtosis (exponential)
t0 = 500
k = 100

skewness = df['Sim_Skewness']
Excess_Kurtosis = df['Sim_Excess_kurtosis']

indices1 = np.arange(t0, len(skewness), k)
indices2 = np.arange(t0, len(Excess_Kurtosis), k)

sim_skewness_piecewise = skewness[indices1]
sim_kurtosis_piecewise = Excess_Kurtosis[indices2]

final_average_sim_skewness = sim_skewness_piecewise.mean()
final_average_sim_kurtosis = sim_kurtosis_piecewise.mean()

print(f"\nSimulated skewness (piecewise_average): {final_average_sim_skewness:.2f}")
print(f"Simulated Excess Kurtosis (piecewise_average): {final_average_sim_kurtosis:.2f}")

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()



