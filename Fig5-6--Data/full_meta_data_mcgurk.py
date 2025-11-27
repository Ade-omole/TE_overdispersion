import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import skew, kurtosis

# Load the data
df = pd.read_excel("path/mcgurk_data.xlsx",
                   sheet_name='Sheet1', header=1)

# Clean up column names
df.columns = df.columns.str.strip()

# Define groups
groups = {
    'DNA transposons': ['POGO', 'P-element', 'M4DM', 'BARI_DM', 'HOBO'],
    'LTR transposons': ['ROO_LTR', 'DM412_LTR', 'NOMAD_LTR', 'MDG1_LTR',
                        'TRANSPAC_LTR', 'BURDOCK_LTR', 'TIRANT_LTR', 'TABOR_LTR'],
    'Non-LTR transposons': ['Jockey', 'FW', 'I_DM', 'DOC', 'DOC6_DM']
}

# Prepare results for each TE individually
results = []
for group_name, te_families in groups.items():
    for te in te_families:
        if te in df.columns:
            values = df[te].dropna().values
            # Filter out copy numbers > 200
            values = values[values <= 200]
            mean_val = values.mean()
            var_val = values.var(ddof=1)  # sample variance
            results.append({
                "Group": group_name,
                "TE": te,
                "Mean": mean_val,
                "Variance": var_val,
                "Variance/Mean": var_val / mean_val if mean_val != 0 else float("inf"),
                "Mean/Variance": mean_val / var_val if var_val != 0 else float("inf"),
                "Skewness": skew(values),
                "Excess Kurtosis": kurtosis(values, fisher=True)
            })
        else:
            print(f"Warning: {te} not found in dataframe columns")

results_df = pd.DataFrame(results)

# Compute group averages of skewness and kurtosis
group_avgs = results_df.groupby("Group")[["Skewness", "Excess Kurtosis"]].mean().reset_index()

# --- Print results ---
print("Per-TE Statistics (Mean, Variance, Ratios, Skewness, Excess Kurtosis):")
print("="*100)
print(results_df.to_string(index=False))


