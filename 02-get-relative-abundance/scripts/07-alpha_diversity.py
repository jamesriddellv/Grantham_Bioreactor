#!/usr/bin/env python
# coding: utf-8

# # This notebook plots alpha diversity statistics based on vOTU GeTMM relative abundance computed from metatranscriptomes

# In[1]:


import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import regex as re
import textwrap
import time
from matplotlib_venn import venn3
import skbio
import matplotlib as mpl


# In[2]:


from print_versions import print_versions
print_versions(globals())


# In[3]:


# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


# In[4]:


# df generated from 06-identify-active-vOTU-and-calc-getmm-abundance.py
df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
df.head()


# In[5]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# # Compare Hill Number package to observed_otus, shannon, and inv_simpson

# In[6]:


df = df.set_index('vOTU')


# In[7]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import observed_otus
import seaborn as sns

# Read data
df_ceil = np.ceil(df)
df_notrans = df


# In[8]:


# Read metadata
metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")

# Simplify treatment labels
metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
metadata = metadata.loc[metadata['treatment'] != 'T']

# Create treatment_timepoint column
treatment_timepoint = metadata['treatment'] + "_" + metadata['timepoint'].astype(str)
df.columns = treatment_timepoint.values

# Make unique column names
def make_unique(columns):
    seen = {}
    new_cols = []
    for col in columns:
        if col not in seen:
            seen[col] = 0
            new_cols.append(col)
        else:
            seen[col] += 1
            new_cols.append(f"{col}.{seen[col]}")
    return new_cols

df_ceil.columns = make_unique(df_ceil.columns)
df_notrans.columns = make_unique(df_notrans.columns)


# In[9]:


# Transpose and convert to numeric matrix
ceil_mtx = df_ceil.T
ceil_mtx = ceil_mtx.apply(pd.to_numeric, errors='coerce')
ceil_mtx


# In[10]:


# Transpose and convert to numeric matrix
df_mtx = df_notrans.T
df_mtx = df_mtx.apply(pd.to_numeric, errors='coerce')
df_mtx


# In[11]:


metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('day')[1]))


# In[12]:


metadata


# In[13]:


from skbio.diversity.alpha import observed_otus, shannon, inv_simpson, hill
import pandas as pd
import numpy as np

# 1. Calculate Classic Alpha Diversity Indices
metadata['ceil_richness'] = ceil_mtx.astype(int).apply(observed_otus, axis=1).values
metadata['ceil_shannon'] = ceil_mtx.astype(int).apply(shannon, axis=1, base=np.e).values 
metadata['ceil_inv_simpson'] = ceil_mtx.astype(int).apply(inv_simpson, axis=1).values

# 2. Calculate Hill Numbers (q=0, 1, 2)
metadata['ceil_hill_q0'] = ceil_mtx.astype(int).apply(lambda x: hill(x, order=0), axis=1).values
metadata['ceil_hill_q1'] = ceil_mtx.astype(int).apply(lambda x: hill(x, order=1), axis=1).values
metadata['ceil_hill_q2'] = ceil_mtx.astype(int).apply(lambda x: hill(x, order=2), axis=1).values

# 3. Print the results for each function
print("--- Alpha Diversity Results (Standard vs. Hill Numbers) ---")
# Create a results dataframe with the correct transformations
results = pd.DataFrame(index=metadata.index)

# q=0: Richness
results['Richness (Observed)'] = metadata['ceil_richness']
results['Hill q=0'] = metadata['ceil_hill_q0']

# q=1: Exponential of Shannon
# We apply np.exp() here so the values actually match the Hill number
results['Exp(Shannon)'] = np.exp(metadata['ceil_shannon'])
results['Hill q=1'] = metadata['ceil_hill_q1']

# q=2: Inverse Simpson
results['Inv Simpson'] = metadata['ceil_inv_simpson']
results['Hill q=2'] = metadata['ceil_hill_q2']

print("--- Results for Each Function ---")
print(results.head())

# Check for equality 
print("\n--- Verification (Equality Check) ---")
print(f"Richness == Hill q0: {np.allclose(results['Richness (Observed)'], results['Hill q=0'])}")
print(f"Exp(Shannon) == Hill q1: {np.allclose(results['Exp(Shannon)'], results['Hill q=1'])}")
print(f"Inv Simpson == Hill q2: {np.allclose(results['Inv Simpson'], results['Hill q=2'])}")


# In[14]:


results


# # Compare hill numbers to rounded-up and non-rounded data to demonstrate artificial inflation from rounding to justify changing the method to reviewers.

# In[15]:


# Compute Hill numbers on non-transformed data
metadata['hill_q0'] = df_mtx.astype(int).apply(lambda x: hill(x, order=0), axis=1).values
metadata['hill_q1'] = df_mtx.astype(int).apply(lambda x: hill(x, order=1), axis=1).values
metadata['hill_q2'] = df_mtx.astype(int).apply(lambda x: hill(x, order=2), axis=1).values


# In[16]:


diff_q0 = metadata['ceil_hill_q0'] - metadata['hill_q0']
diff_q1 = metadata['ceil_hill_q1'] - metadata['hill_q1']
diff_q2 = metadata['ceil_hill_q2'] - metadata['hill_q2']
print("Maximum absolute differences between ceiling method and no transformation using hill orders:")
print(f"Richness of ceiling minus no transformation: {diff_q0.abs().max():.2e}")
print(f"Shannon Index of ceiling minus no transformation:  {diff_q1.abs().max():.2e}")
print(f"Inverse Simpson of ceiling minus no transformation:   {diff_q2.abs().max():.2e}")


# In[17]:


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Define the orders and their corresponding names for titles/filenames
orders = {
    0: "Richness (q=0)",
    1: "Shannon Diversity (q=1)",
    2: "Inverse Simpson (q=2)"
}

# Loop through each order to create separate plots
for q, name in orders.items():
    # 1. Calculate the specific Hill number order
    # Assuming your 'hill' function accepts 'order' as an argument
    metadata['hill_number'] = df_mtx.apply(lambda x: hill(x, order=q), axis=1).values

    # 2. Initialize the plot
    plt.figure(figsize=(6, 6))

    # Points (individual replicates)
    sns.scatterplot(
        data=metadata,
        x='timepoint',
        y='hill_number',
        hue='treatment',
        palette={"C": "#FC9B2D", "U": "#7ACBC3"},
        alpha=0.6 # Slightly lowered alpha to see the mean line better
    )

    # Mean line across timepoints per treatment
    sns.lineplot(
        data=metadata,
        x='timepoint',
        y='hill_number',
        hue='treatment',
        palette={"C": "#FC9B2D", "U": "#7ACBC3"},
        estimator='mean',
        errorbar='ci',
        linewidth=4,
        legend=False
    )

    # 3. Formatting
    plt.title(name, fontsize=16)
    plt.xlabel("Day", fontsize=14)
    plt.ylabel(f"Hill Number (q={q})", fontsize=14)
    
    plt.xticks([0, 7, 14, 21, 35], fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.3)

    # Remove legend and spines
    if plt.gca().get_legend():
        plt.gca().get_legend().remove()
        
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    plt.tight_layout()

    # 4. Save with a dynamic filename
    file_path = f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_hill_order_{q}.pdf'
    plt.savefig(file_path, dpi=300)
    
    print(f"Saved: {file_path}")
    plt.show()


# # Compare side-by-side

# In[ ]:


from skbio.diversity.alpha import observed_otus, shannon, inv_simpson, hill
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define the orders and their corresponding names
orders = {
    0: "Richness (q=0)",
    1: "Shannon Diversity (q=1)",
    2: "Inverse Simpson (q=2)"
}

# Loop through each order to create separate plots showing Ceil vs Raw
for q, name in orders.items():
    
    # 1. Calculate values for both matrices
    # Use copy() to avoid SettingWithCopy warnings
    plot_df = metadata.copy()
    plot_df['Raw'] = df_mtx.apply(lambda x: hill(x, order=q), axis=1).values
    plot_df['Ceiling'] = ceil_mtx.astype(int).apply(lambda x: hill(x, order=q), axis=1).values

    # 2. Melt the dataframe to make it compatible with Seaborn hue/style
    # This turns 'Raw' and 'Ceiling' columns into a single 'Transformation' column
    melted_df = plot_df.melt(
        id_vars=['timepoint', 'treatment'], 
        value_vars=['Raw', 'Ceiling'],
        var_name='Transformation', 
        value_name='hill_value'
    )

    # 3. Initialize the plot
    plt.figure(figsize=(8, 6))

    # Points (individual replicates) - Split by Treatment and Transformation
    sns.scatterplot(
        data=melted_df,
        x='timepoint',
        y='hill_value',
        hue='treatment',
        style='Transformation',
        palette={"C": "#FC9B2D", "U": "#7ACBC3"},
        alpha=0.4,
        legend=True
    )

    # Mean lines - Solid for Raw, Dashed for Ceiling
    sns.lineplot(
        data=melted_df,
        x='timepoint',
        y='hill_value',
        hue='treatment',
        style='Transformation',
        palette={"C": "#FC9B2D", "U": "#7ACBC3"},
        dashes={'Raw': '', 'Ceiling': (2, 2)}, # Solid for Raw, Dashed for Ceil
        estimator='mean',
        errorbar='ci',
        linewidth=3
    )

    # 4. Formatting
    plt.title(f"{name}: Raw vs. Ceiling", fontsize=16)
    plt.xlabel("Day", fontsize=14)
    plt.ylabel(f"Hill Number (q={q})", fontsize=14)
    
    plt.xticks([0, 7, 14, 21, 35], fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, alpha=0.2)

    # Clean up spines
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    # Move legend outside to keep the plot clean
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

    plt.tight_layout()

    # 5. Save
    file_path = f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/Comparison_hill_order_{q}.pdf'
    plt.savefig(file_path, dpi=300, bbox_inches='tight')
    
    print(f"Saved comparison plot for q={q} to: {file_path}")
    plt.show()


# # In addition to alpha diversity, we want to plot and compare the proportion of reads mapping to active vOTUs across treatments and over time.

# In[ ]:


gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt', sep='\t')
gff_df = gff_df[['vOTU', 'gene']]
gff_df.head()

# Read the TSV file without header
counts = pd.read_csv(
    "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED.tsv",
    sep="\t",
    header=None
)

# Assign column names
counts.columns = [
    "gene",
    "STM_0716_E_M_E002",
    "STM_0716_E_M_E003",
    "STM_0716_E_M_E004",
    "STM_0716_E_M_E025",
    "STM_0716_E_M_E027",
    "STM_0716_E_M_E029",
    "STM_0716_E_M_E030",
    "STM_0716_E_M_E031",
    "STM_0716_E_M_E033",
    "STM_0716_E_M_E034",
    "STM_0716_E_M_E035",
    "STM_0716_E_M_E050",
    "STM_0716_E_M_E051",
    "STM_0716_E_M_E052",
    "STM_0716_E_M_E054",
    "STM_0716_E_M_E055",
    "STM_0716_E_M_E056",
    "STM_0716_E_M_E058",
    "STM_0716_E_M_E059",
    "STM_0716_E_M_E060",
    "STM_0716_E_M_E062",
    "STM_0716_E_M_E063",
    "STM_0716_E_M_E064",
    "STM_0716_E_M_E066",
    "STM_0716_E_M_E067",
    "STM_0716_E_M_E068",
    "STM_0716_E_M_E070",
    "STM_0716_E_M_E071",
    "STM_0716_E_M_E072",
    "STM_0716_E_M_E121",
    "STM_0716_E_M_E122",
    "STM_0716_E_M_E123",
    "STM_0716_E_M_E125",
    "STM_0716_E_M_E126",
    "STM_0716_E_M_E127",
    "STM_0716_E_M_E129",
    "STM_0716_E_M_E130",
    "STM_0716_E_M_E131"
]
counts.head()


# In[ ]:


counts = counts.merge(gff_df, on='gene', how='left').dropna()
counts


# In[ ]:


counts_long = counts.melt(id_vars=['gene', 'vOTU'], var_name='Sample', value_name='num_reads_mapped')
counts_long


# In[ ]:


counts_long = counts_long.groupby(['vOTU', 'Sample']).agg({'num_reads_mapped': 'sum'}).reset_index()


# In[ ]:


counts_long = counts_long.merge(replicate_frame, on='Sample', how='left')


# In[ ]:


counts_long = counts_long.loc[counts_long['treatment'] != 'CT']


# In[ ]:


counts_long


# In[ ]:


# keep only active vOTUs
counts_long = counts_long.loc[counts_long['vOTU'].isin(list(df.index))]
counts_long


# In[ ]:


reads_mapped_per_sample = counts_long.groupby(['Sample', 'treatment', 'day', 'replicate']).agg({'num_reads_mapped': 'sum'}).reset_index()
reads_mapped_per_sample


# In[ ]:


metadata.columns


# In[ ]:


reads_mapped_per_sample = reads_mapped_per_sample.merge(metadata[['Sample', 'total reads (R1+R2)']], on='Sample', how='left')


# In[ ]:


reads_mapped_per_sample['prop_mapped'] = reads_mapped_per_sample['num_reads_mapped'] / reads_mapped_per_sample['total reads (R1+R2)'] * 100
reads_mapped_per_sample


# In[ ]:


reads_mapped_per_sample['day'] = reads_mapped_per_sample['day'].astype(int)


# In[ ]:


reads_mapped_per_sample.loc[reads_mapped_per_sample['day'] == 14]


# In[ ]:


metadata = metadata.merge(reads_mapped_per_sample[['Sample', 'prop_mapped']], on='Sample', how='left')
metadata


# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns

# Base plot
fig, ax1 = plt.subplots(figsize=(6, 6))

sns.lineplot(
    data=metadata,
    x='timepoint',
    y='hill_q1',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=1,
    legend=False,
    ax=ax1
)

ax1.set_xlabel("Day", fontsize=14)
ax1.set_ylabel("vOTU Alpha Diversity Exp(H') Index, metaT", fontsize=14)
ax1.set_xticks([0, 7, 14, 21, 35])
ax1.tick_params(axis='both', labelsize=12)
ax1.grid(True, alpha=0.3)

# Remove plot borders
for spine in ax1.spines.values():
    spine.set_visible(False)

# Create second y-axis
ax2 = ax1.twinx()
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='prop_mapped',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    linestyle='--',
    errorbar='ci',
    linewidth=1,
    legend=False,
    ax=ax2
)
ax2.set_ylabel("% of total RNA mapped to active vOTUs", fontsize=14)
ax2.tick_params(axis='y', labelsize=12)

# Final adjustments
fig.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-A_shannon_diversity_with_prop_mapped.pdf', dpi=300)
plt.show()


# # Perform Welch's t-test

# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# ── Aesthetic mappings ────────────────────────────────────────────────────────
color_palette = {"C": "#FC9B2D", "U": "#7ACBC3"}
marker_map = {"0": "o", "7": "s", "14": "^", "21": "s", "35": "P"}

# ── Helper: extract groups (timepoint != 0) ───────────────────────────────────
def get_groups(col):
    cat = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'C')][col].values
    una = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'U')][col].values
    return cat, una

# ── Helper: Welch's one-tailed t-test + Shapiro, prints results ───────────────
def welch_one_tailed(col, h1_greater, label):
    """
    h1_greater : 'C' or 'U'  — which group is hypothesised to be greater
    """
    cat, una = get_groups(col)
    
    # ── log transform for prop_mapped to achieve normal distribution ────────────────────────────────
    if col == 'prop_mapped':
        cat = np.log(cat)
        una = np.log(una)

    # Normality
    wc, pc = stats.shapiro(cat)
    wu, pu = stats.shapiro(una)
    print(f"\n=== Shapiro-Wilk [{label}] ===")
    print(f"  Catechin  W={wc:.4f} p={pc:.4f} {'✓' if pc > 0.05 else '✗'}")
    print(f"  Unamended W={wu:.4f} p={pu:.4f} {'✓' if pu > 0.05 else '✗'}")

    print(f"\n=== Descriptive Stats [{label}] ===")
    print(f"  Catechin  n={len(cat)} mean={np.mean(cat):.4f} SD={np.std(cat, ddof=1):.4f}")
    print(f"  Unamended n={len(una)} mean={np.mean(una):.4f} SD={np.std(una, ddof=1):.4f}")

    # One-tailed Welch's — pass the hypothesised-greater group first
    if h1_greater == 'U':
        t_stat, p_two = stats.ttest_ind(una, cat, equal_var=False)
        h1_label = "Unamended > Catechin"
    else:
        t_stat, p_two = stats.ttest_ind(cat, una, equal_var=False)
        h1_label = "Catechin > Unamended"

    p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2

    print(f"\n=== Welch's t-test (one-tailed: {h1_label}) [{label}] ===")
    print(f"  t = {t_stat:.4f}  p_two = {p_two:.4f}  p_one = {p_one:.4f}")
    sig = p_one < 0.05
    direction = h1_label if sig else f"No evidence that {h1_label}"
    print(f"  → {'Significant' if sig else 'Not significant'}: {direction}")

    return p_one

# ── Run both tests ────────────────────────────────────────────────────────────
p_alpha = welch_one_tailed('hill_q1',   h1_greater='U', label='Alpha Diversity')
p_prop  = welch_one_tailed('prop_mapped', h1_greater='C', label='Prop Mapped')

# ── Helper: draw one stripplot layer ─────────────────────────────────────────
def add_stripplot(ax, plot_df, y_col):
    for tp, m in marker_map.items():
        subset = plot_df[plot_df['timepoint_str'] == tp]
        if subset.empty:
            continue
        if tp == "7":
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker='s', color='black', alpha=0.6, jitter=False, s=9)
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker='x', color='white', alpha=1.0, jitter=False, s=6, linewidth=1)
        else:
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker=m, color='black', alpha=0.6, jitter=True, s=8, linewidth=0.5)

# ── Helper: significance bracket ─────────────────────────────────────────────
def add_bracket(ax, y_vals, p_val):
    y_max = np.nanmax(y_vals)
    y_pos = y_max + y_max * 0.12
    h     = y_max * 0.03
    ax.plot([0, 0, 1, 1], [y_pos, y_pos+h, y_pos+h, y_pos], lw=1.2, c='k')
    ax.text(0.5, y_pos + h + y_max * 0.02, f"p = {p_val:.3e}",
            ha='center', va='bottom', fontsize=9, weight='bold')

# ── Shared plot_df ────────────────────────────────────────────────────────────
plot_df = metadata[~((metadata['treatment'] == 'U') & (metadata['timepoint'] == 0))].copy()
plot_df['timepoint_str'] = plot_df['timepoint'].astype(str)


# In[ ]:


# ── Build 2×2 figure ─────────────────────────────────────────────────────────
sns.set_style("white")

# Create a figure where the first column is 2x the width of the second
fig, axes = plt.subplots(2, 2, figsize=(15, 10), 
                        gridspec_kw={'width_ratios': [2, 0.5]})
(ax_tl, ax_tr), (ax_bl, ax_br) = axes

# ── TOP-LEFT: alpha diversity over time ───────────────────────────────────────
sns.lineplot(data=metadata, x='timepoint', y='hill_q1',
             hue='treatment', palette=color_palette,
             estimator='mean', errorbar='ci', linewidth=1,
             legend=True, ax=ax_tl)

sns.scatterplot(data=metadata, x='timepoint', y='hill_q1',
                hue='treatment', palette=color_palette,
                legend=False, ax=ax_tl, s=50, alpha=0.6)

ax_tl.set_xlabel("Day", fontsize=13)
ax_tl.set_ylabel("vOTU Alpha Diversity Exp(H') Index, metaT", fontsize=12)
ax_tl.set_xticks([0, 7, 14, 21, 35])
# Explicitly show y-axis ticks
ax_tl.tick_params(axis='y', which='both', left=True, labelsize=11) 
ax_tl.grid(True, alpha=0.3)
for sp in ax_tl.spines.values(): sp.set_visible(False)

# ── BOTTOM-LEFT: prop_mapped over time ────────────────────────────────────────
sns.lineplot(data=metadata, x='timepoint', y='prop_mapped',
             hue='treatment', palette=color_palette,
             estimator='mean', errorbar='ci', linewidth=1, linestyle='--',
             legend=True, ax=ax_bl)

sns.scatterplot(data=metadata, x='timepoint', y='prop_mapped',
                hue='treatment', palette=color_palette,
                legend=False, ax=ax_bl, s=50, alpha=0.6)

ax_bl.set_xlabel("Day", fontsize=13)
ax_bl.set_ylabel("% of total RNA mapped to active vOTUs", fontsize=12)
ax_bl.set_xticks([0, 7, 14, 21, 35])
# Explicitly show y-axis ticks
ax_bl.tick_params(axis='y', which='both', left=True, labelsize=11)
ax_bl.grid(True, alpha=0.3)
for sp in ax_bl.spines.values(): sp.set_visible(False)

# ── TOP-RIGHT: boxplot alpha diversity ────────────────────────────────────────
sns.boxplot(data=plot_df, x='treatment', y='hill_q1',
            palette=color_palette, showfliers=False,
            width=0.6, linewidth=1.2, ax=ax_tr)
add_stripplot(ax_tr, plot_df, 'hill_q1')
add_bracket(ax_tr, plot_df['hill_q1'].values, p_alpha)
ax_tr.set_xlabel("Treatment", fontsize=13)
ax_tr.set_ylabel("vOTU Alpha Diversity Exp(H') Index, metaT", fontsize=12)
ax_tr.tick_params(axis='y', left=True) # Ensure ticks show on boxplot

# ── BOTTOM-RIGHT: boxplot prop_mapped ─────────────────────────────────────────
sns.boxplot(data=plot_df, x='treatment', y='prop_mapped',
            palette=color_palette, showfliers=False,
            width=0.6, linewidth=1.2, ax=ax_br)
add_stripplot(ax_br, plot_df, 'prop_mapped')
add_bracket(ax_br, plot_df['prop_mapped'].values, p_prop)
ax_br.set_xlabel("Treatment", fontsize=13)
ax_br.set_ylabel("% of total RNA mapped to active vOTUs", fontsize=12)
ax_br.tick_params(axis='y', left=True) # Ensure ticks show on boxplot

# ── Save ──────────────────────────────────────────────────────────────────────
fig.tight_layout()

plt.show()


# In[ ]:


df_mtx


# ### power analysis

# In[ ]:


# Define functions

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower

def welch_df(s1, s2, n1, n2):
    """Welch-Satterthwaite degrees of freedom."""
    v1, v2 = s1**2 / n1, s2**2 / n2
    return (v1 + v2)**2 / (v1**2 / (n1 - 1) + v2**2 / (n2 - 1))
 
 
def cohens_d(group1, group2):
    """Pooled Cohen's d (uses pooled SD denominator)."""
    n1, n2 = len(group1), len(group2)
    pooled_sd = np.sqrt(
        ((n1 - 1) * group1.std(ddof=1)**2 + (n2 - 1) * group2.std(ddof=1)**2)
        / (n1 + n2 - 2)
    )
    return (group1.mean() - group2.mean()) / pooled_sd if pooled_sd > 0 else np.nan
 
 
def n_for_power(target_power, effect_size, alpha=0.05, alternative='larger'):
    """
    Minimum n per group to achieve target_power for a given effect size.
    Returns np.nan if effect_size is zero or NaN.
    """
    if not np.isfinite(effect_size) or effect_size == 0:
        return np.nan
    analysis = TTestIndPower()
    # solve_power returns float; ceil to next integer
    n = analysis.solve_power(
        effect_size=abs(effect_size),
        alpha=alpha,
        power=target_power,
        alternative=alternative
    )
    return int(np.ceil(n))


# In[ ]:


# ── Part 1: Per-Timepoint Welch Tests ─────────────────────────────────────────

test_results = []
timepoints = sorted(metadata['timepoint'].unique())

for tp in timepoints:
    cat = metadata[(metadata['timepoint'] == tp) & (metadata['treatment'] == 'C')]['prop_mapped'].dropna()
    una = metadata[(metadata['timepoint'] == tp) & (metadata['treatment'] == 'U')]['prop_mapped'].dropna()
    n_cat, n_una = len(cat), len(una)

    p_val, status, df_welch, var_ratio, levene_p, d, power_actual, n_needed_80 = [np.nan]*8

    if n_cat < 2 or n_una < 2:
        status = f"Too few samples (C:{n_cat}, U:{n_una})"
    else:
        s_cat, s_una = cat.std(ddof=1), una.std(ddof=1)
        var_ratio = (s_una**2 / s_cat**2) if s_cat > 0 else np.nan
        _, levene_p = stats.levene(una, cat, center='median')
        
        # Using the helpers (assumed defined: welch_df, cohens_d, n_for_power)
        df_w = welch_df(s_cat, s_una, n_cat, n_una) 
        _, p_v = stats.ttest_ind(cat, una, equal_var=False, alternative='greater')
        eff_d = cohens_d(cat, una) 
        
        if np.isfinite(eff_d) and np.isfinite(df_w):
            t_crit = stats.t.ppf(0.95, df=df_w)
            ncp = abs(eff_d) / np.sqrt(1/n_cat + 1/n_una)
            power_actual = 1 - stats.nct.cdf(t_crit, df=df_w, nc=ncp)
        
        n_80 = n_for_power(0.80, eff_d, alternative='larger')
        
        # Map back to variables for dict
        p_val, df_welch, d, n_needed_80 = p_v, df_w, eff_d, n_80

    test_results.append({
        'Timepoint': str(tp), 'n_C': n_cat, 'n_U': n_una, 'var_ratio': var_ratio,
        'levene_p': levene_p, 'df_welch': df_welch, 'cohens_d': d,
        'power_actual': power_actual, 'n_needed_80': n_needed_80, 'p_raw': p_val, 'Notes': status
    })

results_tp = pd.DataFrame(test_results)

# ── Part 2: Combination & FDR Correction ──────────────────────────────────────

results_combined = results_tp # no more combination

# Correct for multiple testing across all P-values (Per-TP and GLS Aggregate)
valid_mask = results_combined['p_raw'].notna()
if valid_mask.any():
    reject, adj_p, _, _ = multipletests(results_combined.loc[valid_mask, 'p_raw'], method='fdr_bh')
    results_combined['p_adj'] = np.nan
    results_combined.loc[valid_mask, 'p_adj'] = adj_p
    results_combined['Significant'] = False
    results_combined.loc[valid_mask, 'Significant'] = reject

# Formatting and Display
pd.options.display.float_format = '{:.3f}'.format
print("\n=== FINAL RESULTS (Per-TP Welch) ===")
cols_to_show = ['Timepoint', 'n_C', 'n_U', 'cohens_d', 'power_actual', 'n_needed_80', 'p_raw', 'p_adj', 'Significant']
print(results_combined[cols_to_show])


# In[ ]:


results_combined.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/power_analysis_prop_mapped_table.csv', index=False)


# In[ ]:


# ── Part 1: Per-Timepoint Welch Tests ─────────────────────────────────────────

test_results = []
timepoints = sorted(metadata['timepoint'].unique())

for tp in timepoints:
    cat = metadata[(metadata['timepoint'] == tp) & (metadata['treatment'] == 'C')]['hill_q1'].dropna()
    una = metadata[(metadata['timepoint'] == tp) & (metadata['treatment'] == 'U')]['hill_q1'].dropna()
    n_cat, n_una = len(cat), len(una)

    p_val, status, df_welch, var_ratio, levene_p, d, power_actual, n_needed_80 = [np.nan]*8

    if n_cat < 2 or n_una < 2:
        status = f"Too few samples (C:{n_cat}, U:{n_una})"
    else:
        s_cat, s_una = cat.std(ddof=1), una.std(ddof=1)
        var_ratio = (s_una**2 / s_cat**2) if s_cat > 0 else np.nan
        _, levene_p = stats.levene(una, cat, center='median')
        
        # Using the helpers (assumed defined: welch_df, cohens_d, n_for_power)
        df_w = welch_df(s_cat, s_una, n_cat, n_una) 
        _, p_v = stats.ttest_ind(una, cat, equal_var=False, alternative='greater')
        eff_d = cohens_d(una, cat) 
        
        if np.isfinite(eff_d) and np.isfinite(df_w):
            t_crit = stats.t.ppf(0.95, df=df_w)
            ncp = abs(eff_d) / np.sqrt(1/n_cat + 1/n_una)
            power_actual = 1 - stats.nct.cdf(t_crit, df=df_w, nc=ncp)
        
        n_80 = n_for_power(0.80, eff_d, alternative='larger')
        
        # Map back to variables for dict
        p_val, df_welch, d, n_needed_80 = p_v, df_w, eff_d, n_80

    test_results.append({
        'Timepoint': str(tp), 'n_C': n_cat, 'n_U': n_una, 'var_ratio': var_ratio,
        'levene_p': levene_p, 'df_welch': df_welch, 'cohens_d': d,
        'power_actual': power_actual, 'n_needed_80': n_needed_80, 'p_raw': p_val, 'Notes': status
    })

results_tp = pd.DataFrame(test_results)

# ── Part 2: Combination & FDR Correction ──────────────────────────────────────

results_combined = results_tp # no more combination

# Correct for multiple testing across all P-values (Per-TP and GLS Aggregate)
valid_mask = results_combined['p_raw'].notna()
if valid_mask.any():
    reject, adj_p, _, _ = multipletests(results_combined.loc[valid_mask, 'p_raw'], method='fdr_bh')
    results_combined['p_adj'] = np.nan
    results_combined.loc[valid_mask, 'p_adj'] = adj_p
    results_combined['Significant'] = False
    results_combined.loc[valid_mask, 'Significant'] = reject

# Formatting and Display
pd.options.display.float_format = '{:.3f}'.format
print("\n=== FINAL RESULTS (Per-TP Welch) ===")
cols_to_show = ['Timepoint', 'n_C', 'n_U', 'cohens_d', 'power_actual', 'n_needed_80', 'p_raw', 'p_adj', 'Significant']
print(results_combined[cols_to_show])


# In[ ]:


results_combined.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/power_analysis_hill_q1_table.csv', index=False)


# # Test time dependency

# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np

# Assuming 'metadata' is your DataFrame
# 1. Statistical Testing (Pooling Treatments)

# Spearman Correlation (Trend over time)
rho, p_spearman = stats.spearmanr(metadata['timepoint'], metadata['hill_q1'])

# Kruskal-Wallis (Difference between any timepoints)
# Grouping data by timepoint for the test
groups = [group['hill_q1'].values for name, group in metadata.groupby('timepoint')]
h_stat, p_kruskal = stats.kruskal(*groups)

print(f"=== Statistical Results (Pooled Treatments) ===")
print(f"Spearman Correlation: rho = {rho:.4f}, p = {p_spearman:.4f}")
print(f"Kruskal-Wallis Test: H = {h_stat:.4f}, p = {p_kruskal:.4f}")

# 2. Generate the Boxplot
plt.figure(figsize=(8, 6))

# Use a neutral color for combined treatments
sns.boxplot(data=metadata, x='timepoint', y='hill_q1', color='lightgrey', showfliers=False)
sns.stripplot(data=metadata, x='timepoint', y='hill_q1', color='black', alpha=0.5, jitter=True)

# Add title and labels
plt.title('Alpha Diversity (Hill q1) Over Time (Treatments Combined)', fontsize=14)
plt.xlabel('Timepoint (Days)', fontsize=12)
plt.ylabel('Alpha Diversity (Hill q1)', fontsize=12)

# Annotate with p-value from Kruskal-Wallis
plt.text(0.5, 0.95, f'Kruskal-Wallis p = {p_kruskal:.4f}\nSpearman rho = {rho:.4f} (p = {p_spearman:.4f})', 
         horizontalalignment='center', verticalalignment='top', 
         transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.5))

plt.tight_layout()
plt.savefig('alpha_diversity_over_time.png', dpi=300)


# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd

# ── Aesthetic mappings ────────────────────────────────────────────────────────
color_palette = {"C": "#FC9B2D", "U": "#7ACBC3"}
marker_map = {"0": "o", "7": "s", "14": "^", "21": "s", "35": "P"}

# ── 1. TIME NORMALIZATION (Z-SCORING) ─────────────────────────────────────────
# We calculate Z-scores within each timepoint to "remove" the day-to-day variance.
# For prop_mapped, we log-transform first (as in your original logic) then Z-score.

metadata['hill_q1_z'] = metadata.groupby('timepoint')['hill_q1'].transform(
    lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
)

metadata['prop_mapped_log'] = np.log(metadata['prop_mapped'])
metadata['prop_mapped_z'] = metadata.groupby('timepoint')['prop_mapped_log'].transform(
    lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
)

# ── Verification: Kruskal-Wallis test on Normalized Data ─────────────────────

# 1. Test Alpha Diversity Z-scores across time
groups_alpha_z = [group['hill_q1_z'].values for name, group in metadata.groupby('timepoint')]
h_alpha, p_alpha_time = stats.kruskal(*groups_alpha_z)

# 2. Test Prop Mapped Z-scores across time
groups_prop_z = [group['prop_mapped_z'].values for name, group in metadata.groupby('timepoint')]
h_prop, p_prop_time = stats.kruskal(*groups_prop_z)

print("=== Time Dependency Check (Post-Normalization) ===")
print(f"Alpha Diversity Z-score vs Time: H = {h_alpha:.4f}, p = {p_alpha_time:.4f}")
print(f"Prop Mapped Z-score vs Time:     H = {h_prop:.4f}, p = {p_prop_time:.4f}")

if p_alpha_time > 0.05 and p_prop_time > 0.05:
    print("\nSUCCESS: Time dependency has been statistically removed.")
else:
    print("\nWARNING: Some variance remains; check for timepoints with very low sample sizes.")

# ── Helper: extract groups (timepoint != 0) ───────────────────────────────────
def get_groups(col):
    cat = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'C')][col].values
    una = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'U')][col].values
    return cat, una

# ── Helper: Welch's one-tailed t-test ─────────────────────────────────────────
def welch_one_tailed(col, h1_greater, label):
    cat, una = get_groups(col)
    
    # Normality on normalized data
    wc, pc = stats.shapiro(cat)
    wu, pu = stats.shapiro(una)
    print(f"\n=== Shapiro-Wilk (Normalized) [{label}] ===")
    print(f"  Catechin  W={wc:.4f} p={pc:.4f} {'✓' if pc > 0.05 else '✗'}")
    print(f"  Unamended W={wu:.4f} p={pu:.4f} {'✓' if pu > 0.05 else '✗'}")

    # One-tailed Welch's
    if h1_greater == 'U':
        t_stat, p_two = stats.ttest_ind(una, cat, equal_var=False)
        h1_label = "Unamended > Catechin"
    else:
        t_stat, p_two = stats.ttest_ind(cat, una, equal_var=False)
        h1_label = "Catechin > Unamended"

    p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2
    print(f"\n=== Welch's t-test on Z-Scores (one-tailed: {h1_label}) [{label}] ===")
    print(f"  t = {t_stat:.4f}  p_one = {p_one:.4f}")
    
    return p_one

# ── Run tests on Normalized Data ──────────────────────────────────────────────
p_alpha = welch_one_tailed('hill_q1_z',   h1_greater='U', label='Alpha Diversity')
p_prop  = welch_one_tailed('prop_mapped_z', h1_greater='C', label='Prop Mapped')

# ── Helper: draw one stripplot layer ─────────────────────────────────────────
def add_stripplot(ax, plot_df, y_col):
    for tp, m in marker_map.items():
        subset = plot_df[plot_df['timepoint_str'] == tp]
        if subset.empty: continue
        if tp == "7":
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker='s', color='black', alpha=0.6, jitter=False, s=9)
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker='x', color='white', alpha=1.0, jitter=False, s=6, linewidth=1)
        else:
            sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                          marker=m, color='black', alpha=0.6, jitter=True, s=8, linewidth=0.5)

# ── Helper: significance bracket ─────────────────────────────────────────────
def add_bracket(ax, y_vals, p_val):
    y_max = np.nanmax(y_vals)
    y_pos = y_max + 0.5 # Small offset for Z-score scale
    h     = 0.1
    ax.plot([0, 0, 1, 1], [y_pos, y_pos+h, y_pos+h, y_pos], lw=1.2, c='k')
    ax.text(0.5, y_pos + h + 0.1, f"p = {p_val:.4f}",
            ha='center', va='bottom', fontsize=9, weight='bold')

# ── Update plot_df with Z-scored columns ──────────────────────────────────────
plot_df = metadata[~((metadata['treatment'] == 'U') & (metadata['timepoint'] == 0))].copy()
plot_df['timepoint_str'] = plot_df['timepoint'].astype(str)

# ── Build 2×2 figure ─────────────────────────────────────────────────────────
sns.set_style("white")
fig, axes = plt.subplots(2, 2, figsize=(12, 7), gridspec_kw={'width_ratios': [2, 0.5]})
(ax_tl, ax_tr), (ax_bl, ax_br) = axes

# ── LEFT PANELS: Raw Data Over Time (To show the actual trend) ────────────────
# Top-Left: Raw Alpha
sns.lineplot(data=metadata, x='timepoint', y='hill_q1', hue='treatment', palette=color_palette, ax=ax_tl)
sns.scatterplot(data=metadata, x='timepoint', y='hill_q1', hue='treatment', palette=color_palette, ax=ax_tl, s=50, alpha=0.6, legend=False)
ax_tl.set_ylabel("Raw Alpha Diversity", fontsize=12)

# Bottom-Left: Raw Prop Mapped
sns.lineplot(data=metadata, x='timepoint', y='prop_mapped', hue='treatment', palette=color_palette, ax=ax_bl, linestyle='--')
sns.scatterplot(data=metadata, x='timepoint', y='prop_mapped', hue='treatment', palette=color_palette, ax=ax_bl, s=50, alpha=0.6, legend=False)
ax_bl.set_ylabel("Raw % RNA Mapped", fontsize=12)

# ── RIGHT PANELS: Normalized Data (Showing Treatment effect ONLY) ─────────────
# Top-Right: Z-scored Alpha
sns.boxplot(data=plot_df, x='treatment', y='hill_q1_z', palette=color_palette, showfliers=False, width=0.6, ax=ax_tr)
add_stripplot(ax_tr, plot_df, 'hill_q1_z')
add_bracket(ax_tr, plot_df['hill_q1_z'].values, p_alpha)
ax_tr.axhline(0, color='black', lw=0.8, ls='--', alpha=0.4) # Mean of each timepoint
ax_tr.set_ylabel("Alpha Diversity (Z-score by Day)", fontsize=12)

# Bottom-Right: Z-scored Prop Mapped
sns.boxplot(data=plot_df, x='treatment', y='prop_mapped_z', palette=color_palette, showfliers=False, width=0.6, ax=ax_br)
add_stripplot(ax_br, plot_df, 'prop_mapped_z')
add_bracket(ax_br, plot_df['prop_mapped_z'].values, p_prop)
ax_br.axhline(0, color='black', lw=0.8, ls='--', alpha=0.4) # Mean of each timepoint
ax_br.set_ylabel("Prop Mapped (Z-score by Day)", fontsize=12)

# ── Clean up spines and save ──────────────────────────────────────────────────
for ax in [ax_tl, ax_bl]:
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig(
    '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-A_combined_panel.pdf',
    dpi=300, bbox_inches='tight'
)
plt.show()


# In[ ]:


plot_df.to_csv('data/2A_diversity_indices.csv', index=False)


# # Post-hoc power analysis

# In[ ]:


import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.power import TTestIndPower

# --- Updated Helper Functions ---

def cohens_d(group1, group2):
    """Calculates Cohen's d for two independent groups."""
    n1, n2 = len(group1), len(group2)
    # Pooled SD
    pooled_sd = np.sqrt(
        ((n1 - 1) * np.var(group1, ddof=1) + (n2 - 1) * np.var(group2, ddof=1)) 
        / (n1 + n2 - 2)
    )
    return (np.mean(group1) - np.mean(group2)) / pooled_sd if pooled_sd > 0 else 0

def post_hoc_power(group1, group2, alpha=0.05, alternative='larger'):
    """
    Calculates the achieved power (post-hoc) based on actual effect size 
    and sample sizes observed in the data.
    """
    d = abs(cohens_d(group1, group2))
    n1, n2 = len(group1), len(group2)
    ratio = n2 / n1
    
    if d == 0:
        return 0.0
        
    analysis = TTestIndPower()
    power = analysis.solve_power(
        effect_size=d,
        nobs1=n1,
        ratio=ratio,
        alpha=alpha,
        alternative=alternative
    )
    return power

# --- Updated Analysis Wrapper ---

def analyze_with_power(col, h1_greater, label):
    cat, una = get_groups(col)
    n_cat, n_una = len(cat), len(una)
    
    # 1. Normality Check
    wc, pc = stats.shapiro(cat)
    wu, pu = stats.shapiro(una)
    
    # 2. t-test Logic
    # If h1_greater is 'U', we test if Unamended > Catechin
    group_a = una if h1_greater == 'U' else cat
    group_b = cat if h1_greater == 'U' else una
    h1_label = f"{'Unamended' if h1_greater == 'U' else 'Catechin'} > Others"

    t_stat, p_two = stats.ttest_ind(group_a, group_b, equal_var=False)
    p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2
    
    # 3. Post-Hoc Power & Effect Size
    d = cohens_d(group_a, group_b)
    achieved_power = post_hoc_power(group_a, group_b, alternative='larger')
    
    print(f"\n=== Results for {label} (Z-Scores) ===")
    print(f"  Shapiro p: Cat={pc:.4f}, Una={pu:.4f}")
    print(f"  Hypothesis: {h1_label}")
    print(f"  t-stat: {t_stat:.4f} | p-value: {p_one:.4f}")
    print(f"  --- Post-Hoc Power Analysis ---")
    print(f"  Observed Cohen's d: {d:.4f}")
    print(f"  Sample sizes: n_cat={n_cat}, n_una={n_una}")
    print(f"  Achieved Power: {achieved_power:.4f}")
    
    return p_one

# --- Execution ---
p_alpha = analyze_with_power('hill_q1_z', h1_greater='U', label='Alpha Diversity')
p_prop  = analyze_with_power('prop_mapped_z', h1_greater='C', label='Prop Mapped')


# In[ ]:


import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.power import TTestIndPower

def calculate_metrics(df, col, h1_greater, timepoint):
    """Calculates t-test and power for a specific timepoint subset."""
    # Subset data
    subset = df[df['timepoint'] == timepoint]
    cat = subset[subset['treatment'] == 'C'][col].values
    una = subset[subset['treatment'] == 'U'][col].values
    
    n_cat, n_una = len(cat), len(una)
    
    # Handle cases with insufficient data
    if n_cat < 2 or n_una < 2:
        return {
            'Timepoint': timepoint, 'n_cat': n_cat, 'n_una': n_una,
            't_stat': np.nan, 'p_val': np.nan, 'Cohen_d': np.nan, 'Power': np.nan
        }

    # Define groups based on hypothesis
    group_a = una if h1_greater == 'U' else cat
    group_b = cat if h1_greater == 'U' else una

    # Statistics
    t_stat, p_two = stats.ttest_ind(group_a, group_b, equal_var=False)
    p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2
    
    # Cohen's d
    pooled_sd = np.sqrt(((n_cat-1)*cat.std(ddof=1)**2 + (n_una-1)*una.std(ddof=1)**2) / (n_cat+n_una-2))
    d = (group_a.mean() - group_b.mean()) / pooled_sd if pooled_sd > 0 else 0
    
    # Post-hoc Power
    analysis = TTestIndPower()
    pw = analysis.solve_power(effect_size=abs(d), nobs1=len(group_a), 
                              ratio=len(group_b)/len(group_a), alpha=0.05, alternative='larger')

    return {
        'Timepoint': timepoint,
        'n_cat': n_cat,
        'n_una': n_una,
        't_stat': round(t_stat, 3),
        'p_val': round(p_one, 4),
        'Cohen_d': round(d, 3),
        'Power': round(pw, 4)
    }

# --- Run Analysis for All Timepoints ---

results_alpha = []
results_prop = []
timepoints = sorted([t for t in metadata['timepoint'].unique() if t != 0])

for tp in timepoints:
    results_alpha.append(calculate_metrics(metadata, 'hill_q1', 'U', tp))
    results_prop.append(calculate_metrics(metadata, 'prop_mapped_log', 'C', tp))

# Create Tables
df_alpha = pd.DataFrame(results_alpha)
df_prop = pd.DataFrame(results_prop)

print("\n--- POST-HOC POWER: ALPHA DIVERSITY (H1: U > C) ---")
print(df_alpha.to_string(index=False))

print("\n--- POST-HOC POWER: PROP MAPPED (H1: C > U---")
print(df_prop.to_string(index=False))


# In[ ]:


df_alpha


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 003-alpha-diversity-sensitivity-analysis.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts/07-alpha_diversity')


# In[ ]:




