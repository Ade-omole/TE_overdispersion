import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

# Set publication-ready style
plt.style.use('default')
sns.set_palette("Set2")

# --- Load and process data (your existing code) ---
df = pd.read_excel("path/mcgurk_data.xlsx",
                   sheet_name='Sheet1', header=1)
df.columns = df.columns.str.strip()

groups = {
    'DNA transposons': ['POGO', 'P-element', 'M4DM', 'BARI_DM', 'HOBO'],
    'LTR transposons': ['ROO_LTR', 'DM412_LTR', 'NOMAD_LTR', 'MDG1_LTR',
                        'TRANSPAC_LTR', 'BURDOCK_LTR', 'TIRANT_LTR', 'TABOR_LTR'],
    'Non-LTR transposons': ['Jockey', 'FW', 'I_DM', 'DOC', 'DOC6_DM']
}

# Calculate statistics with copy number filtering
results = []
for group_name, te_families in groups.items():
    for te in te_families:
        if te in df.columns:
            values = df[te].dropna().values
            # Filter out copy numbers > 200
            values = values[values <= 200]
            if len(values) > 1:
                mean_val = values.mean()
                var_val = values.var(ddof=1)
                dispersion_index = var_val / mean_val if mean_val > 0 else np.nan
                results.append({
                    "Group": group_name,
                    "TE": te,
                    "Mean": mean_val,
                    "Variance": var_val,
                    "Dispersion_Index": dispersion_index,
                    "Log_Mean": np.log10(mean_val) if mean_val > 0 else np.nan,
                    "Log_Variance": np.log10(var_val) if var_val > 0 else np.nan,
                    "N_samples": len(values)
                })

results_df = pd.DataFrame(results)
results_df = results_df.dropna()

# Define colors for groups
colors = {'DNA transposons': '#2E8B57',
          'LTR transposons': '#FF6B35',
          'Non-LTR transposons': '#4472C4'}

# Create figure with subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# ========================================
# PLOT A: Variance vs Mean Scatter Plot (ENHANCED)
# ========================================
for group in groups.keys():
    subset = results_df[results_df["Group"] == group]
    ax1.scatter(subset["Mean"], subset["Variance"],
                label=group, color=colors[group], s=200, alpha=0.9,  # Increased size and alpha
                edgecolors='black', linewidth=3)  # Thicker, darker edge

# Poisson line and overdispersion zones
min_val = results_df["Mean"].min() * 0.5
max_val = results_df["Mean"].max() * 2
x_line = np.logspace(np.log10(min_val), np.log10(max_val), 100)

ax1.plot(x_line, x_line, "k--", label="Poisson (Variance = Mean)", linewidth=3.5, alpha=0.9)
ax1.fill_between(x_line, x_line, x_line * 10, alpha=0.12, color='red',
                 label='Moderate overdispersion (2-10×)')

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Mean Copy Number (log10)", fontsize=24, fontweight='bold')
ax1.set_ylabel("Variance in Copy Number (log10)", fontsize=24, fontweight='bold')

# Enhanced legend with increased font size and boldness
legend1 = ax1.legend(frameon=True, fancybox=True, shadow=True,
                    loc='upper left', framealpha=0.95, edgecolor='black',
                    prop={'size': 14, 'weight': 'bold'})
legend1.get_frame().set_linewidth(1.2)

ax1.grid(alpha=0.4, which="both", linestyle=":", linewidth=0.8)

# Enhanced tick parameters with larger, bolder numbers
ax1.tick_params(axis='both', which='major', labelsize=20,
                width=2, length=6)
ax1.tick_params(axis='both', which='minor', labelsize=20,
                width=1.5, length=4)

# Make axis tick labels bold
for label in ax1.get_xticklabels():
    label.set_fontweight('bold')
for label in ax1.get_yticklabels():
    label.set_fontweight('bold')

# ========================================
# PLOT B: Individual TE Elements Ranked
# ========================================
all_elements = results_df.sort_values('Dispersion_Index')
colors_list = [colors[group] for group in all_elements['Group']]

bars = ax2.bar(range(len(all_elements)), all_elements['Dispersion_Index'],
               color=colors_list, alpha=0.8, edgecolor='white', linewidth=1)

# Add reference lines
ax2.axhline(y=1, color='black', linestyle='--', linewidth=3, alpha=0.9,
            label='Poisson expectation')
ax2.axhline(y=10, color='red', linestyle='--', linewidth=3, alpha=0.9,
            label='High overdispersion threshold')

ax2.set_yscale('log')
ax2.set_xlabel('Transposable elements', fontsize=20, fontweight='bold')
ax2.set_ylabel('Variance/Mean (log10)', fontsize=22, fontweight='bold')

# Add selective annotations for extreme values and every few elements
annotation_indices = []
# Add highest and lowest few elements
annotation_indices.extend(range(3))  # First 3 (lowest)
annotation_indices.extend(range(len(all_elements) - 3, len(all_elements)))  # Last 3 (highest)
# Add elements with very high dispersion
for i, (idx, row) in enumerate(all_elements.iterrows()):
    if row['Dispersion_Index'] > 1:
        annotation_indices.append(i)

# Remove duplicates and sort
annotation_indices = sorted(list(set(annotation_indices)))

# ENHANCED ELEMENT NAME ANNOTATIONS - LARGER AND BOLDER
for i in annotation_indices:
    row = all_elements.iloc[i]
    ax2.annotate(row['TE'], (i, row['Dispersion_Index']),
                 xytext=(0, 8), textcoords='offset points',
                 fontsize=14, rotation=90, ha='left', fontweight='bold')  # Increased from 10 to 14

# Create perfect legend for TE classes
legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[group], alpha=0.8,
                                 edgecolor='black', linewidth=1, label=group)
                   for group in groups.keys()]

# Add reference lines to legend
legend_elements.extend([
    plt.Line2D([0], [0], color='black', linestyle='--', linewidth=4,
               label='Poisson expectation'),
    plt.Line2D([0], [0], color='red', linestyle='--', linewidth=4,
               label='High overdispersion threshold')
])

# Enhanced legend with increased font size and boldness
legend2 = ax2.legend(handles=legend_elements, loc='upper left',
                    bbox_to_anchor=(0, 0.92),  # Move it slightly down from top
                    frameon=True, fancybox=True, shadow=True,
                    framealpha=0.95, edgecolor='black',
                    prop={'size': 14, 'weight': 'bold'})
legend2.get_frame().set_linewidth(1.2)

ax2.grid(alpha=0.4, axis='y', linestyle=':', linewidth=0.8)

# Enhanced tick parameters with larger, bolder numbers
ax2.tick_params(axis='both', which='major', labelsize=18,
                width=2, length=6)
ax2.tick_params(axis='both', which='minor', labelsize=18,
                width=1.5, length=4)

# Make axis tick labels bold
for label in ax2.get_xticklabels():
    label.set_fontweight('bold')
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')

ax2.set_xticks([])  # Remove x-axis ticks as they're not meaningful here

plt.tight_layout()
plt.show()