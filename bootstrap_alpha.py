import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys




# ================= Configuration Area =================
INPUT_FILE = sys.argv[1]  # Please replace with your filename
N_BOOTSTRAP = 2000            # Number of bootstrap resamples, recommended 1000-5000
RANDOM_SEED = 42              # Random seed to ensure reproducibility
TARGET_LEN_COL = 'length_of_ctcf'
NEUTRAL_LEN_COL = 'length_of_gene'
# ======================================================




def calculate_alpha(df):
    """
    Calculate Alpha value
    Logic: Normalize length for each row first, then aggregate counts across all rows, finally calculate alpha
    """
    # 1. Length normalization factor
    # Scale neutral region counts to "how many mutations would there be if the neutral region were only as long as the target"
    # Factor = L_target / L_neutral
    scale_factor = df[TARGET_LEN_COL] / df[NEUTRAL_LEN_COL]
    
    # 2. Calculate scaled neutral region counts (keep as float for precision)
    # D_neu_scaled = D_neu_raw * (L_tar / L_neu)
    d_neu_scaled = df['dsyn'] * scale_factor
    p_neu_scaled = df['psyn'] * scale_factor
    
    # 3. Aggregate counts across all sites (Summation is crucial for MK tests)
    sum_d_tar = df['dctcf'].sum()
    sum_p_tar = df['pctcf'].sum()
    sum_d_neu = d_neu_scaled.sum()
    sum_p_neu = p_neu_scaled.sum()
    
    # 4. Prevent division by zero
    if sum_p_neu == 0 or sum_d_neu == 0:
        return np.nan
    
    # 5. Calculate Alpha
    # Alpha = 1 - (P_tar / P_neu_scaled) / (D_tar / D_neu_scaled)
    poly_ratio = sum_p_tar / sum_p_neu
    div_ratio = sum_d_tar / sum_d_neu
    
    if div_ratio == 0:
        return np.nan
        
    alpha = 1 - (poly_ratio / div_ratio)
    return alpha




# ================= Main Program =================




print(f"Reading data file: {INPUT_FILE} ...")
try:
    df = pd.read_csv(INPUT_FILE, sep='\t')
except FileNotFoundError:
    print(f"Error: File {INPUT_FILE} not found")
    exit(1)




# Data cleaning: Remove rows with NaN in required columns
required_cols = ['dctcf', 'pctcf', 'dsyn', 'psyn', TARGET_LEN_COL, NEUTRAL_LEN_COL]
df = df.dropna(subset=required_cols)




# Filter out rows with length 0 to prevent division by zero
df = df[(df[TARGET_LEN_COL] > 0) & (df[NEUTRAL_LEN_COL] > 0)]




if len(df) == 0:
    print("Error: No data remaining after cleaning.")
    exit(1)




print(f"Number of valid sites: {len(df)}")




# 1. Calculate Observed Value (Observed Alpha)
alpha_obs = calculate_alpha(df)
print(f"\nObserved Alpha (Observed α): {alpha_obs:.4f}")




# 2. Bootstrap Resampling
print(f"Performing Bootstrap resampling ({N_BOOTSTRAP} times)...")
np.random.seed(RANDOM_SEED)
bootstrap_alphas = []




for i in range(N_BOOTSTRAP):
    # Sampling with replacement (frac=1.0 means sample size equals original data size)
    sample_df = df.sample(n=len(df), replace=True)
    alpha_star = calculate_alpha(sample_df)
    bootstrap_alphas.append(alpha_star)




# Clean NaN values (if any)
bootstrap_alphas = [x for x in bootstrap_alphas if not np.isnan(x)]




if len(bootstrap_alphas) < 100:
    print("Warning: Too few valid Bootstrap samples; results may be unreliable.")
else:
    # 3. Calculate 95% Confidence Interval
    ci_lower = np.percentile(bootstrap_alphas, 2.5)
    ci_upper = np.percentile(bootstrap_alphas, 97.5)
    
    print(f"\n=== Analysis Results ===")
    print(f"Observed α: {alpha_obs:.4f}")
    print(f"95% CI:     [{ci_lower:.4f}, {ci_upper:.4f}]")
    
   
    # 4. Plotting
    plt.figure(figsize=(10, 6))
    plt.hist(bootstrap_alphas, bins=40, color='skyblue', edgecolor='black', alpha=0.7, density=False)
    
    # Mark observed value
    plt.axvline(alpha_obs, color='red', linestyle='--', linewidth=2, label=f'Observed α={alpha_obs:.3f}')
    # Mark CI
    plt.axvline(ci_lower, color='green', linestyle=':', linewidth=2, label=f'95% CI Lower={ci_lower:.3f}')
    plt.axvline(ci_upper, color='green', linestyle=':', linewidth=2, label=f'95% CI Upper={ci_upper:.3f}')
    # Mark zero line
    plt.axvline(0, color='black', linestyle='-', linewidth=1)
    
    plt.xlabel('Alpha Estimate (Proportion of Adaptation)')
    plt.ylabel('Frequency')
    plt.title(f'Bootstrap Distribution of Alpha (N={len(bootstrap_alphas)})')
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    output_img = sys.argv[2]
    plt.savefig(output_img, dpi=300)
    print(f"\nHistogram saved as: {output_img}")