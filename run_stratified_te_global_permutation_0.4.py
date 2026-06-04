#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import pybedtools
import sys
import os


# ================= 1. CONFIGURATION =================


# [EDIT 1] Input CBS file (Must have the header as shown)
INPUT_CBS_FILE = sys.argv[1] 


# [EDIT 2] Superpool File 
# CRITICAL: This file MUST also have the SAME 4 annotation columns (GC_Bin, etc.)
# It can have a header too, but the column names must match exactly.
SUPERPOOL_FILE = "superpool_25bp_example.bed" 


# [EDIT 3] TE Annotation File (Standard BED or RepeatMasker format)
TE_FILE = "rmsk.sort.example.bed" 


# [EDIT 4] Output Prefix
OUTPUT_PREFIX = INPUT_CBS_FILE.split('merged_robust')[0] + "_0.2_strat_result"


# [EDIT 5] Stratification Column Names (Must match your file header exactly)
STRATA_COLS = ['GC_Bin', 'Promoter_Bin', 'Mappability_Bin', 'Distance_Bin']


# Parameters
MOTIF_LEN = 25
OVERLAP_THRESH = 0.20
MIN_OVERLAP_BP = int(MOTIF_LEN * OVERLAP_THRESH)
N_PERM = 10000  # 10,000 iterations for robust P-values
SEED = 42


np.random.seed(SEED)


print(f"🚀 Starting Strictly Stratified TE Permutation Test...")
print(f"   Input CBS: {INPUT_CBS_FILE}")
print(f"   Superpool: {SUPERPOOL_FILE}")
print(f"   Stratification Keys: {STRATA_COLS}")
print(f"   Overlap Threshold: >= {MIN_OVERLAP_BP}bp ({OVERLAP_THRESH*100}%)")
print(f"   Permutations: {N_PERM}")


# ================= 2. LOAD DATA WITH HEADERS =================


print("\n📊 Loading data with headers...")


try:
    # Load CBS with header
    cbs_df = pd.read_csv(INPUT_CBS_FILE, sep='\t', header=0)
    
    # Load Superpool with header
    sp_df = pd.read_csv(SUPERPOOL_FILE, sep='\t', header=0)
    
    te = pybedtools.BedTool(TE_FILE)
except Exception as e:
    print(f"❌ File loading failed: {e}")
    print("   Please ensure both input files have the correct header row.")
    sys.exit(1)


# Validate Columns
missing_cbs_cols = [col for col in STRATA_COLS if col not in cbs_df.columns]
missing_sp_cols = [col for col in STRATA_COLS if col not in sp_df.columns]


if missing_cbs_cols:
    print(f"❌ Error: CBS file missing columns: {missing_cbs_cols}")
    sys.exit(1)
if missing_sp_cols:
    print(f"❌ Error: Superpool file missing columns: {missing_sp_cols}")
    print("   Your Superpool file MUST be annotated with the same 4 columns as the CBS file.")
    sys.exit(1)


print(f"   ✅ Columns verified. Found {len(cbs_df):,} CBS regions.")


# ================= 3. CREATE STRATA KEYS =================


# Create a unique key by combining the 4 columns (e.g., "high_intergenic_low_near")
# Fill NaNs with 'NA' just in case to avoid errors during string conversion
cbs_df[STRATA_COLS] = cbs_df[STRATA_COLS].fillna('NA')
sp_df[STRATA_COLS] = sp_df[STRATA_COLS].fillna('NA')


cbs_df['strata_key'] = cbs_df[STRATA_COLS].astype(str).agg('_'.join, axis=1)
sp_df['strata_key'] = sp_df[STRATA_COLS].astype(str).agg('_'.join, axis=1)


print(f"   Found {cbs_df['strata_key'].nunique()} unique strata in CBS.")
print(f"   Found {sp_df['strata_key'].nunique()} unique strata in Superpool.")


# Prepare BedTools objects (only need chrom, start, end for geometry)
cbs_bed = pybedtools.BedTool.from_dataframe(cbs_df[['Chr', 'Start', 'End']])
sp_bed = pybedtools.BedTool.from_dataframe(sp_df[['Chr', 'Start', 'End']])


# ================= 4. PRE-COMPUTE TE STATUS FOR SUPERPOOL =================


print(f"\n⚡ Pre-computing TE overlap for Superpool (One-time cost)...")


sp_intersect = sp_bed.intersect(te, wo=True)


# Base dataframe for merging results
sp_base = sp_df[['Chr', 'Start', 'End']].copy()
sp_base['coord_key'] = sp_base['Chr'].astype(str) + ":" + sp_base['Start'].astype(str) + "-" + sp_base['End'].astype(str)


if sp_intersect.count() == 0:
    sp_base['overlap_len'] = 0
else:
    df_int = sp_intersect.to_dataframe(names=None)
    overlap_col_idx = df_int.shape[1] - 1
    
    # Group by coords and sum overlap
    grouped = df_int.iloc[:, overlap_col_idx].groupby([df_int.iloc[:, 0], df_int.iloc[:, 1], df_int.iloc[:, 2]]).sum()
    grouped_df = grouped.reset_index()
    grouped_df.columns = ['Chr', 'Start', 'End', 'overlap_len']
    grouped_df['coord_key'] = grouped_df['Chr'].astype(str) + ":" + grouped_df['Start'].astype(str) + "-" + grouped_df['End'].astype(str)
    
    merged = sp_base.merge(grouped_df[['coord_key', 'overlap_len']], on='coord_key', how='left')
    sp_base['overlap_len'] = merged['overlap_len'].fillna(0).astype(int)


# Apply Threshold
sp_base['is_te_assoc'] = (sp_base['overlap_len'] >= MIN_OVERLAP_BP).astype(int)


# Attach status back to full superpool dataframe (keeping strata_key)
sp_df['is_te_assoc'] = sp_base['is_te_assoc'].values


# Group Superpool by strata for fast sampling
# Result: dict where key=strata_name, value=list of 0/1 status
sp_groups = sp_df.groupby('strata_key')['is_te_assoc'].apply(list).to_dict()
sp_counts = {k: len(v) for k, v in sp_groups.items()}


print(f"   Superpool annotated. Global TE fraction: {sp_df['is_te_assoc'].mean():.4f}")


# ================= 5. CALCULATE OBSERVED STATISTIC =================


print(f"\n🧮 Calculating Observed TE overlap for Real CBSs...")


cbs_intersect = cbs_bed.intersect(te, wo=True)


if cbs_intersect.count() == 0:
    cbs_df['is_te_assoc'] = 0
else:
    df_cbs_int = cbs_intersect.to_dataframe(names=None)
    overlap_col_idx = df_cbs_int.shape[1] - 1
    
    grouped_cbs = df_cbs_int.iloc[:, overlap_col_idx].groupby([df_cbs_int.iloc[:, 0], df_cbs_int.iloc[:, 1], df_cbs_int.iloc[:, 2]]).sum()
    
    cbs_base = cbs_df[['Chr', 'Start', 'End']].copy()
    cbs_base['coord_key'] = cbs_base['Chr'].astype(str) + ":" + cbs_base['Start'].astype(str) + "-" + cbs_base['End'].astype(str)
    
    grouped_df = grouped_cbs.reset_index()
    grouped_df.columns = ['Chr', 'Start', 'End', 'overlap_len']
    grouped_df['coord_key'] = grouped_df['Chr'].astype(str) + ":" + grouped_df['Start'].astype(str) + "-" + grouped_df['End'].astype(str)
    
    merged_cbs = cbs_base.merge(grouped_df[['coord_key', 'overlap_len']], on='coord_key', how='left')
    merged_cbs['overlap_len'] = merged_cbs['overlap_len'].fillna(0).astype(int)
    
    cbs_df['is_te_assoc'] = (merged_cbs['overlap_len'] >= MIN_OVERLAP_BP).astype(int)


obs_total_hits = cbs_df['is_te_assoc'].sum()
n_total_cbs = len(cbs_df)
obs_fraction = obs_total_hits / n_total_cbs if n_total_cbs > 0 else 0


print(f"   Observed TE-associated Count: {obs_total_hits}")
print(f"   Observed Fraction: {obs_fraction:.4f}")


# ================= 6. RUN STRATIFIED PERMUTATION TEST =================


# Count CBSs per stratum
cbs_group_counts = cbs_df.groupby('strata_key').size().to_dict()


# Validation: Ensure every stratum in CBS has enough candidates in Superpool
valid_strata = []
for k, n_needed in cbs_group_counts.items():
    if k not in sp_counts:
        print(f"   ⚠️ Warning: Stratum '{k}' found in CBS but NOT in Superpool. Skipping this stratum.")
        continue
    if sp_counts[k] < n_needed:
        print(f"   ❌ CRITICAL ERROR: Stratum '{k}' has {n_needed} CBSs but only {sp_counts[k]} candidates in Superpool.")
        print(f"   Cannot perform valid permutation test. Please expand your Superpool.")
        sys.exit(1)
    valid_strata.append(k)


if not valid_strata:
    print("   ❌ No valid strata found for permutation. Exiting.")
    sys.exit(1)


print(f"   Running {N_PERM} stratified permutations across {len(valid_strata)} valid strata...")


perm_fractions = np.zeros(N_PERM)


for i in range(N_PERM):
    total_hits_perm = 0
    total_n_perm = 0
    
    for strat in valid_strata:
        n_sample = cbs_group_counts[strat]
        pool_array = np.array(sp_groups[strat])
        
        # Stratified Sampling: Sample ONLY from the matching stratum
        sampled = np.random.choice(pool_array, size=n_sample, replace=False)
        total_hits_perm += sampled.sum()
        total_n_perm += n_sample
    
    perm_fractions[i] = total_hits_perm / total_n_perm if total_n_perm > 0 else 0
    
    if (i + 1) % 2000 == 0:
        print(f"      Progress: {i+1}/{N_PERM}")


# Calculate P-value
# One-sided test: Is Observed significantly HIGHER than random?
p_val = (np.sum(perm_fractions >= obs_fraction) + 1) / (N_PERM + 1)


print(f"\n📈 Final Stratified Results:")
print(f"   - Observed Fraction: {obs_fraction:.6f}")
print(f"   - Random Mean (Expected): {np.mean(perm_fractions):.6f}")
print(f"   - P-Value: {p_val:.6f}")
print(f"   - N_CBS: {n_total_cbs}")


# ================= 7. SAVE OUTPUT =================


output_file = f"{OUTPUT_PREFIX}_stats.txt"


with open(output_file, 'w') as f:
    f.write("Metric\tValue\n")
    f.write(f"Observed_Fraction\t{obs_fraction}\n")
    f.write(f"Superpool_Mean\t{np.mean(perm_fractions)}\n")
    f.write(f"P_Value\t{p_val}\n")
    f.write(f"N_CBS\t{n_total_cbs}\n")
    f.write(f"Observed_Count\t{obs_total_hits}\n")
    f.write(f"Stratification_Cols\t{','.join(STRATA_COLS)}\n")


print(f"✅ Results saved to: {output_file}")
