import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import skew, kurtosis

# Load data
df = pd.read_excel("path/mcgurk_data.xlsx",
                   sheet_name='Sheet1', header=1)
df.columns = df.columns.str.strip()

groups = {
    'DNA transposons': ['POGO', 'P-element', 'M4DM', 'BARI_DM', 'HOBO'],
    'LTR transposons': ['ROO_LTR', 'DM412_LTR', 'NOMAD_LTR', 'MDG1_LTR',
                        'TRANSPAC_LTR', 'BURDOCK_LTR', 'TIRANT_LTR', 'TABOR_LTR'],
    'Non-LTR transposons': ['Jockey', 'FW', 'I_DM', 'DOC', 'DOC6_DM']
}

# Define colors for groups
colors = {'DNA transposons': '#2E8B57',
          'LTR transposons': '#FF6B35',
          'Non-LTR transposons': '#4472C4'}

# Get all TE elements
all_elements = []
element_groups = []
for group_name, te_list in groups.items():
    for te in te_list:
        if te in df.columns:
            all_elements.append(te)
            element_groups.append(group_name)

# Split into two groups - adjusted split
n_elements = len(all_elements)
# Move one more element to second plot for better balance
split_point = (n_elements // 2) - 1

# PLOT 1: First portion of elements
elements_1 = all_elements[:split_point]
groups_1 = element_groups[:split_point]

cols = 4
rows = int(np.ceil(len(elements_1) / cols))

fig, axes = plt.subplots(rows, cols, figsize=(16, 4 * rows))
if rows == 1:
    axes = axes.reshape(1, -1)
axes = axes.flatten()

for i, (element, group) in enumerate(zip(elements_1, groups_1)):
    ax = axes[i]

    data = df[element].dropna()
    mean_val = data.mean()
    var_val = data.var(ddof=1)
    skew_val = skew(data)
    kurt_val = kurtosis(data, fisher=True)

    sns.histplot(data, bins=15, kde=True, color=colors[group],
                 edgecolor='black', linewidth=1, alpha=0.7,
                 line_kws={'linewidth': 2, 'color': 'darkred'}, ax=ax)

    ax.set_title(f'{element}', fontsize=14, fontweight='bold')

    # Only show xlabel for bottom row
    if i >= (rows - 1) * cols or i >= len(elements_1) - cols:
        ax.set_xlabel('Copy Number', fontsize=12, fontweight='bold')
    else:
        ax.set_xlabel('')

    # Only show ylabel for leftmost column
    if i % cols == 0:
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('')

    ax.tick_params(axis='both', labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    textstr = '\n'.join((
        f'μ = {mean_val:.1f}',
        f'σ² = {var_val:.1f}',
        f'Skew = {skew_val:.2f}',
        f'Ex. Kurt = {kurt_val:.2f}',
    ))
    ax.text(0.98, 0.95, textstr, transform=ax.transAxes,
            fontsize=9, fontweight='bold',
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

# Hide empty subplots
for i in range(len(elements_1), len(axes)):
    axes[i].set_visible(False)

plt.tight_layout()
plt.show()

# PLOT 2: Second portion of elements
elements_2 = all_elements[split_point:]
groups_2 = element_groups[split_point:]

rows = int(np.ceil(len(elements_2) / cols))

fig, axes = plt.subplots(rows, cols, figsize=(16, 4 * rows))
if rows == 1:
    axes = axes.reshape(1, -1)
axes = axes.flatten()

for i, (element, group) in enumerate(zip(elements_2, groups_2)):
    ax = axes[i]

    data = df[element].dropna()
    mean_val = data.mean()
    var_val = data.var(ddof=1)
    skew_val = skew(data)
    kurt_val = kurtosis(data, fisher=True)

    sns.histplot(data, bins=15, kde=True, color=colors[group],
                 edgecolor='black', linewidth=1, alpha=0.7,
                 line_kws={'linewidth': 2, 'color': 'darkred'}, ax=ax)

    ax.set_title(f'{element}', fontsize=14, fontweight='bold')

    # Only show xlabel for bottom row
    if i >= (rows - 1) * cols or i >= len(elements_2) - cols:
        ax.set_xlabel('Copy Number', fontsize=12, fontweight='bold')
    else:
        ax.set_xlabel('')

    # Only show ylabel for leftmost column
    if i % cols == 0:
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('')

    ax.tick_params(axis='both', labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    textstr = '\n'.join((
        f'μ = {mean_val:.1f}',
        f'σ² = {var_val:.1f}',
        f'Skew = {skew_val:.2f}',
        f'Ex. Kurt = {kurt_val:.2f}',
    ))
    ax.text(0.98, 0.95, textstr, transform=ax.transAxes,
            fontsize=9, fontweight='bold',
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

# Hide empty subplots
for i in range(len(elements_2), len(axes)):
    axes[i].set_visible(False)

# Remove legend for second plot
# Create legend for groups
# fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.95),
#            fontsize=12, frameon=True, fancybox=True, shadow=True,
#            framealpha=0.95, edgecolor='black', prop={'weight': 'bold'})

plt.tight_layout()
plt.show()

# PLOT 3: Special comparison of P-element, DOC6_DM, and MDG1_LTR
special_elements = ['P-element', 'DOC6_DM', 'MDG1_LTR']
special_groups = ['DNA transposons', 'Non-LTR transposons', 'LTR transposons']

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, (element, group) in enumerate(zip(special_elements, special_groups)):
    ax = axes[i]

    data = df[element].dropna()
    mean_val = data.mean()
    var_val = data.var(ddof=1)
    skew_val = skew(data)
    kurt_val = kurtosis(data, fisher=True)

    sns.histplot(data, bins=15, kde=True, color=colors[group],
                 edgecolor='black', linewidth=1, alpha=0.7,
                 line_kws={'linewidth': 2, 'color': 'darkred'}, ax=ax)

    ax.set_title(f'{element}', fontsize=16, fontweight='bold')
    ax.set_xlabel('Copy Number', fontsize=14, fontweight='bold')

    if i == 0:
        ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    else:
        ax.set_ylabel('')

    ax.tick_params(axis='both', labelsize=12)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    textstr = '\n'.join((
        f'μ = {mean_val:.1f}',
        f'σ² = {var_val:.1f}',
        f'Skew = {skew_val:.2f}',
        f'Ex. Kurt = {kurt_val:.2f}',
    ))
    ax.text(0.98, 0.95, textstr, transform=ax.transAxes,
            fontsize=11, fontweight='bold',
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

# Print summary statistics for all elements
print("\nSummary Statistics for All TE Elements:")
print("=" * 80)
print(f"{'Element':<15} {'Group':<20} {'Mean':<8} {'Var':<8} {'Skew':<8} {'Kurt':<8}")
print("-" * 80)

for element, group in zip(all_elements, element_groups):
    data = df[element].dropna()
    mean_val = data.mean()
    var_val = data.var(ddof=1)
    skew_val = skew(data)
    kurt_val = kurtosis(data, fisher=True)

    print(f"{element:<15} {group:<20} {mean_val:>6.1f} {var_val:>6.1f} {skew_val:>6.2f} {kurt_val:>6.2f}")