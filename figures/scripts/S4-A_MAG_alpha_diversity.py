#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Plot MAG alpha diversity from expression data in McGivern et al. 2025


# In[2]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import matplotlib as mpl
# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


from print_versions import print_versions
print_versions(globals())


# In[3]:


metaT = pd.read_excel('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_abundance_table.xlsx')
metaT.head()


# In[4]:


# Drop Condensed Tannin samples
metaT = metaT.drop(columns=['STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                        'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                        'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                        'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127'])


# In[5]:


# drop zeroes
metaT = metaT[metaT.iloc[:, 1:-1].sum(axis=1) > 0]

# Reshape to long format
long_df = metaT.melt(id_vars=["MAG", "GTDB"], 
                  var_name="Sample", 
                  value_name="abundance")
long_df


# In[6]:


# add replicates
replicate_frame = pd.DataFrame({
    'STM_0716_E_M_E002_unamended_0': 1,
    'STM_0716_E_M_E003_unamended_0': 2,
       'STM_0716_E_M_E004_unamended_0': 3,
    'STM_0716_E_M_E025_unamended_14': 1,
       'STM_0716_E_M_E027_unamended_14': 2,
    'STM_0716_E_M_E029_CT_14': 1,
       'STM_0716_E_M_E030_CT_14': 2,
    'STM_0716_E_M_E031_CT_14': 3,
       'STM_0716_E_M_E033_catechin_14': 1,
    'STM_0716_E_M_E034_catechin_14': 2,
       'STM_0716_E_M_E035_catechin_14': 3,
    'STM_0716_E_M_E050_unamended_21': 1,
       'STM_0716_E_M_E051_unamended_21': 2,
    'STM_0716_E_M_E052_unamended_21': 3,
       'STM_0716_E_M_E054_CT_21': 1,
    'STM_0716_E_M_E055_CT_21': 2,
       'STM_0716_E_M_E056_CT_21': 3,
    'STM_0716_E_M_E058_catechin_21': 1,
       'STM_0716_E_M_E059_catechin_21': 2,
    'STM_0716_E_M_E060_catechin_21': 3,
       'STM_0716_E_M_E062_unamended_35': 1,
    'STM_0716_E_M_E063_unamended_35': 2,
       'STM_0716_E_M_E064_unamended_35': 3,
    'STM_0716_E_M_E066_CT_35': 1,
       'STM_0716_E_M_E067_CT_35': 2,
    'STM_0716_E_M_E068_CT_35': 3,
       'STM_0716_E_M_E070_catechin_35': 1,
    'STM_0716_E_M_E071_catechin_35': 2,
       'STM_0716_E_M_E072_catechin_35': 3,
    'STM_0716_E_M_E121_unamended_7': 1,
       'STM_0716_E_M_E122_unamended_7': 2,
    'STM_0716_E_M_E123_unamended_7': 3,
       'STM_0716_E_M_E125_CT_7': 1,
    'STM_0716_E_M_E126_CT_7': 2,
       'STM_0716_E_M_E127_CT_7': 3,
    'STM_0716_E_M_E129_catechin_7': 1,
       'STM_0716_E_M_E130_catechin_7': 2,
    'STM_0716_E_M_E131_catechin_7': 3
}.items())

replicate_frame['Sample'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[0])
replicate_frame['treatment'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[1])
replicate_frame['day'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[2])
replicate_frame.columns = ['sample_treatment_day', 'replicate', 'Sample', 'treatment', 'day']
replicate_frame['treatment_day_replicate'] = replicate_frame['treatment'] + '_' + replicate_frame['day'].astype(str) + '_' + replicate_frame['replicate'].astype(str)


# In[7]:


melted_df = long_df.merge(replicate_frame, on='Sample', how='left')
melted_df['replicate'] = melted_df['replicate'].astype(str)
melted_df.head()


# In[8]:


melted_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/getmms_REV_active_MAGs_melted.tsv', sep='\t', index=False)


# In[9]:


melted_df.MAG.nunique()


# In[10]:


active_MAGs = melted_df


# In[11]:


active_MAGs['day'] = active_MAGs['day'].astype(int)


# In[12]:


active_MAGs


# In[13]:


active_MAGs[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = active_MAGs['GTDB'].str.split(';', expand=True)


# In[14]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

active_MAGs['highest_host_tax_rank'] = active_MAGs.apply(get_highest_resolution, axis=1)


# In[15]:


active_MAGs.loc[active_MAGs['highest_host_tax_rank'] == 'g__Methanosarcina'].MAG.unique()


# In[16]:


min_value = active_MAGs['abundance'].min()
shift_value = 1 - min_value
active_MAGs['log2_abundance'] = active_MAGs['abundance'].apply(lambda x: np.log2(x+shift_value))


# In[17]:


active_MAGs.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_abundance_with_metadata.tsv', sep='\t', index=False)


# In[18]:


active_MAGs.head()


# In[19]:


data_mtx = active_MAGs.pivot_table(values='abundance', index='Sample', columns='MAG')


# In[20]:


# Read metadata
metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")

# Simplify treatment labels
metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('y')[1]))


# In[21]:


metadata = metadata.loc[metadata['treatment'] != 'T']


# In[22]:


# Plot MAG diversity statistics
from skbio.diversity.alpha import hill
metadata['hill_q1'] = data_mtx.apply(lambda x: hill(x, order=1), axis=1).values


# In[23]:


metadata


# In[24]:


# Plot individual points and a mean line per treatment
plt.figure(figsize=(6, 6))

# Points (individual replicates)
sns.scatterplot(
    data=metadata,
    x='timepoint',
    y='hill_q1',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    alpha=1
)

# Mean line across timepoints per treatment
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='hill_q1',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=2,
    legend=False  # Prevent duplicate legends
)

# plt.title("Shannon Diversity")
plt.xlabel("Day", fontsize=14)
plt.ylabel("Hill Number (order=1)", fontsize=14)
treatment_colors = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, alpha=0.3)

# Remove legend
plt.legend().remove()

plt.tight_layout()
# plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_MAG_alpha_diversity.pdf', dpi=300)
plt.show()



# In[30]:


# include significance test
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

# ── Run test ────────────────────────────────────────────────────────────
p_alpha = welch_one_tailed('hill_q1',   h1_greater='U', label='Alpha Diversity')

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


# In[31]:


# ── Build 1×2 figure ─────────────────────────────────────────────────────────
sns.set_style("white")

# Create a figure with 1 row and 2 columns
# Width ratio 2:1 gives the lineplot more space for the time axis
fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 5), 
                                        gridspec_kw={'width_ratios': [2, 1]})

# ── LEFT: Alpha diversity over time (Lineplot) ───────────────────────────────
sns.lineplot(data=metadata, x='timepoint', y='hill_q1',
             hue='treatment', palette=color_palette,
             estimator='mean', errorbar='ci', linewidth=1.5,
             legend=True, ax=ax_left)

sns.scatterplot(data=metadata, x='timepoint', y='hill_q1',
                hue='treatment', palette=color_palette,
                legend=False, ax=ax_left, s=60, alpha=0.6)

ax_left.set_xlabel("Day", fontsize=13)
ax_left.set_ylabel("vOTU Alpha Diversity Exp(H'), metaT", fontsize=12)
ax_left.set_xticks([0, 7, 14, 21, 35])
ax_left.tick_params(axis='y', which='both', left=True, labelsize=11) 
ax_left.grid(True, alpha=0.3)

# Despine for a cleaner look
for sp in ax_left.spines.values(): sp.set_visible(False)

# ── RIGHT: Alpha diversity comparison (Boxplot) ──────────────────────────────
sns.boxplot(data=plot_df, x='treatment', y='hill_q1',
            palette=color_palette, showfliers=False,
            width=0.6, linewidth=1.2, ax=ax_right)

# Adding points and stats (assuming your custom helper functions exist)
add_stripplot(ax_right, plot_df, 'hill_q1')
add_bracket(ax_right, plot_df['hill_q1'].values, p_alpha)

ax_right.set_xlabel("Treatment", fontsize=13)
ax_right.set_ylabel("Alpha Diversity Index", fontsize=12)
ax_right.tick_params(axis='y', left=True)

# ── Save ──────────────────────────────────────────────────────────────────────
fig.tight_layout()
plt.savefig(
    '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4-A_MAG_alpha_diversity.pdf',
    dpi=300, bbox_inches='tight'
)
plt.show()


# In[33]:


import seaborn as sns
import matplotlib.pyplot as plt

# ── 1. TIME NORMALIZATION (Z-SCORING) ─────────────────────────────────────────
# We calculate Z-scores within each timepoint to "remove" the day-to-day variance.
metadata['hill_q1_z'] = metadata.groupby('timepoint')['hill_q1'].transform(
    lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
)
p_alpha = welch_one_tailed('hill_q1_z',   h1_greater='U', label='Alpha Diversity')

# NOTE: Make sure plot_df is generated/filtered AFTER the z-score step above
# so that it also contains the 'hill_q1_z' column for the right-hand plot.

# ── Build 1×2 figure ─────────────────────────────────────────────────────────
sns.set_style("white")

# Create a figure with 1 row and 2 columns
# Width ratio 2:1 gives the lineplot more space for the time axis
fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 5), 
                                        gridspec_kw={'width_ratios': [2, 1]})

# ── LEFT: Alpha diversity over time (Original Values) ────────────────────────
sns.lineplot(data=metadata, x='timepoint', y='hill_q1',
             hue='treatment', palette=color_palette,
             estimator='mean', errorbar='ci', linewidth=1.5,
             legend=True, ax=ax_left)

sns.scatterplot(data=metadata, x='timepoint', y='hill_q1',
                hue='treatment', palette=color_palette,
                legend=False, ax=ax_left, s=60, alpha=0.6)

ax_left.set_xlabel("Day", fontsize=13)
ax_left.set_ylabel("vOTU Alpha Diversity Exp(H'), metaT", fontsize=12)
ax_left.set_xticks([0, 7, 14, 21, 35])
ax_left.tick_params(axis='y', which='both', left=True, labelsize=11) 
ax_left.grid(True, alpha=0.3)

# Despine for a cleaner look
for sp in ax_left.spines.values(): sp.set_visible(False)

# ── RIGHT: Alpha diversity comparison (Z-scored) ─────────────────────────────
sns.boxplot(data=plot_df, x='treatment', y='hill_q1_z',
            palette=color_palette, showfliers=False,
            width=0.6, linewidth=1.2, ax=ax_right)

# Adding points and stats (assuming your custom helper functions exist)
add_stripplot(ax_right, plot_df, 'hill_q1_z')
add_bracket(ax_right, plot_df['hill_q1_z'].values, p_alpha)

ax_right.set_xlabel("Treatment", fontsize=13)
ax_right.set_ylabel("Alpha Diversity Index (Z-score)", fontsize=12)
ax_right.tick_params(axis='y', left=True)

# ── Save ──────────────────────────────────────────────────────────────────────
fig.tight_layout()
plt.savefig(
    '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4-A_MAG_alpha_diversity_split_axis.pdf',
    dpi=300, bbox_inches='tight'
)
plt.show()


# In[1]:


get_ipython().system('jupyter nbconvert --to script 05-MAG_abundance_hill_number.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S2-MAG_alpha_diversity')


# In[ ]:




