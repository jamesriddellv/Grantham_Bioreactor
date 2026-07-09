#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


# In[2]:


df = pd.read_excel('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/supplementary_data_6_wraf108.xlsx', sheet_name='catechin_degradation_genes')
df.head()


# In[3]:


completeness = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/CheckM2_MAGs/quality_report.tsv', usecols=['Name', 'Completeness', 'Contamination'], sep='\t')
completeness = completeness.rename(columns={'Name': 'MAG'})
completeness.head()


# In[4]:


df = df.merge(completeness, on='MAG', how='left')
df.head()


# In[5]:


df = df.loc[(df['Completeness'] >50) & (df['Contamination'] < 10)]
len(df)


# In[6]:


df[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = df['GTDB'].str.split(';', expand=True)

# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return 'no_assignment' # If no valid value is found, return NaN

df['highest_host_tax_rank'] = df.apply(get_highest_resolution, axis=1)


# In[7]:


len(df)


# In[8]:


with open('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_list_no_prefix.txt', 'r') as f:
    active_MAGs = [i.strip() for i in f.readlines()]
active_MAGs[:5]


# In[9]:


# filter for only active MAGs
df = df.loc[df['MAG'].isin(active_MAGs)]


# In[10]:


df.MAG.nunique()


# In[11]:


len(df)


# In[12]:


df['enzyme'].unique()


# # STITCH 2026-01-26: Figure out what additional polyphenol metabolisms were occurring

# In[13]:


# big file, 3.4Gb. Load only CAMPER annotations
all_annotations = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/reactorEMERGE_annotations.tsv', sep='\t', usecols=['Unnamed: 0', 'fasta', 'scaffold', 'gene_position', 'camper_hits', 'camper_rank', 'camper_bitScore', 'camper_id', 'camper_definition', 'camper_search_type'])
all_annotations.shape


# In[14]:


all_annotations.fasta.tolist()[:5]


# In[15]:


active_camper_annotations = all_annotations.loc[(all_annotations['fasta'].isin(df.MAG.unique())) & ~(all_annotations['camper_hits'].isna())]
active_camper_annotations.rename(columns={'Unnamed: 0': 'gene'},inplace=True)


# In[16]:


# Figure out which of these are UPSTREAM of phloroglucinol degradation. Check the polyphenol CAMPER map.
active_camper_annotations.camper_definition.value_counts()


# In[17]:


len(active_camper_annotations)


# # Now figure out if these pathways were active, and if they were upregulated in catechin-amended samples.

# In[18]:


# load per-gene activity for vOTUs, across all samples, isolate these genes.
genes_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/geTMMs_ge5_1X_REV_26samples.txt', sep='\t')
genes_df


# In[19]:


camper_genes_df = genes_df.loc[genes_df['gene'].isin(active_camper_annotations['gene'].unique())]
camper_genes_df


# In[20]:


#### filtered out are those that were inactive.


# In[21]:


camper_genes_df = camper_genes_df.merge(active_camper_annotations, on='gene', how='left')


# In[22]:


camper_genes_df['MAG'] = camper_genes_df['fasta']


# In[23]:


camper_genes_df


# # Split MAGs into having catechin, and non-catechin-degrading phloroglucinol degraders

# In[24]:


camper_genes_df['camper_hits'].unique()


# In[25]:


# Get all MAGs that have either PGR or pgthAB (excluding those with FCR, PHY, CHI)
catechin_degraders = set(df.loc[df['enzyme'].isin(['FCR', 'PHY', 'CHI'])]['MAG'].unique())
pgr_mags = set(df.loc[df['enzyme'] == 'PGR']['MAG'].unique()) - catechin_degraders
pgthab_mags = set(df.loc[df['enzyme'] == 'pgthAB']['MAG'].unique()) - catechin_degraders

# Get all unique MAGs (union of both)
all_mags = pgr_mags | pgthab_mags


# In[26]:


len(catechin_degraders), len(pgr_mags), len(pgthab_mags), len(all_mags)


# In[27]:


PGR_only_df = df.loc[(df['MAG'].isin(pgr_mags)) & (df['enzyme'] == 'PGR')].sort_values(by=['MAG'])
PGR_total_per_sample = PGR_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
PGR_total_per_sample.columns = ['Sample', 'total_getmm']

replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
PGR_total_per_sample = PGR_total_per_sample.merge(replicate_frame, on='Sample', how='left')
PGR_total_per_sample.head()


# In[28]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Calculate mean and 95% confidence interval for each treatment-day combination
summary = PGR_total_per_sample.groupby(['treatment', 'day'])['total_getmm'].agg([
    ('mean', 'mean'),
    ('std', 'std'),
    ('count', 'count')
]).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Define colors for treatments
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Plot each treatment
for treatment in summary['treatment'].unique():
    treatment_data = summary[summary['treatment'] == treatment]
    
    # Plot replicate points
    replicate_data = PGR_total_per_sample[PGR_total_per_sample['treatment'] == treatment]
    ax.scatter(replicate_data['day'], replicate_data['total_getmm'],
              color=colors.get(treatment, 'gray'), alpha=0.6, s=50, zorder=3)
    
    # Plot mean line (no markers)
    ax.plot(treatment_data['day'], treatment_data['mean'], 
            linewidth=2.5,
            label=treatment, color=colors.get(treatment, 'gray'))
    
# Formatting
ax.set_xlabel('Day', fontsize=12, fontweight='bold')
ax.set_ylabel('Total getmm', fontsize=12, fontweight='bold')
ax.set_title('PGR Total Over Time by Treatment', 
             fontsize=14, fontweight='bold')
ax.legend(title='Treatment', fontsize=10)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xticks([0, 7, 14, 21, 35])

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10A-PGR_total_over_time_non-catechin-degraders.pdf', dpi=300)
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(summary.to_string(index=False))


# In[29]:


from scipy import stats

# 1. Filter the data for Days 14-35 and specific treatments
mask_days = PGR_total_per_sample['day'].isin([14, 21, 35])

catechin_vals = PGR_total_per_sample[
    (mask_days) & (PGR_total_per_sample['treatment'] == 'catechin')
]['total_getmm']

unamended_vals = PGR_total_per_sample[
    (mask_days) & (PGR_total_per_sample['treatment'] == 'unamended')
]['total_getmm']

# 2. Perform Welch's t-test (equal_var=False)
# 'greater' specifies the one-sided alternative: catechin > unamended
t_stat, p_val = stats.ttest_ind(catechin_vals, 
                                unamended_vals, 
                                equal_var=False, 
                                alternative='greater')

# 3. Report results
print(f"--- Welch's t-test Results (Days 14-35) ---")
print(f"Alternative Hypothesis: Catechin > Unamended")
print(f"t-statistic: {t_stat:.4f}")
print(f"p-value:     {p_val:.4e}")

if p_val < 0.05:
    print("Result: Statistically Significant (p < 0.05)")
else:
    print("Result: Not Statistically Significant")


# In[30]:


pgth_only_df = df.loc[(df['MAG'].isin(pgthab_mags)) & (df['enzyme'] == 'pgthAB')].sort_values(by=['MAG'])
print(pgth_only_df.MAG.nunique())
pgth_total_per_sample = pgth_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
pgth_total_per_sample.columns = ['Sample', 'total_getmm']
pgth_total_per_sample = pgth_total_per_sample.merge(replicate_frame, on='Sample', how='left')


# Calculate mean and 95% confidence interval for each treatment-day combination
summary = pgth_total_per_sample.groupby(['treatment', 'day'])['total_getmm'].agg([
    ('mean', 'mean'),
    ('std', 'std'),
    ('count', 'count')
]).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Define colors for treatments
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Plot each treatment
for treatment in summary['treatment'].unique():
    treatment_data = summary[summary['treatment'] == treatment]
    
    # Plot replicate points
    replicate_data = pgth_total_per_sample[pgth_total_per_sample['treatment'] == treatment]
    ax.scatter(replicate_data['day'], replicate_data['total_getmm'],
              color=colors.get(treatment, 'gray'), alpha=0.6, s=50, zorder=3)
    
    # Plot mean line (no markers)
    ax.plot(treatment_data['day'], treatment_data['mean'], 
            linewidth=2.5,
            label=treatment, color=colors.get(treatment, 'gray'))
    
# Formatting
ax.set_xlabel('Day', fontsize=12, fontweight='bold')
ax.set_ylabel('Total getmm', fontsize=12, fontweight='bold')
ax.set_title('pgth Total Over Time by Treatment', 
             fontsize=14, fontweight='bold')
ax.legend(title='Treatment', fontsize=10)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xticks([0, 7, 14, 21, 35])

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10B-pgth_total_over_time_non-catechin-degraders.pdf', dpi=300)
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(summary.to_string(index=False))


# In[31]:


from scipy import stats

# 1. Filter the data for Days 14-35 and specific treatments
mask_days = pgth_total_per_sample['day'].isin([14, 21, 35])

catechin_vals = pgth_total_per_sample[
    (mask_days) & (pgth_total_per_sample['treatment'] == 'catechin')
]['total_getmm']

unamended_vals = pgth_total_per_sample[
    (mask_days) & (pgth_total_per_sample['treatment'] == 'unamended')
]['total_getmm']

# 2. Perform Welch's t-test (equal_var=False)
# 'greater' specifies the one-sided alternative: catechin > unamended
t_stat, p_val = stats.ttest_ind(catechin_vals, 
                                unamended_vals, 
                                equal_var=False, 
                                alternative='greater')

# 3. Report results
print(f"--- Welch's t-test Results (Days 14-35) ---")
print(f"Alternative Hypothesis: Catechin > Unamended")
print(f"t-statistic: {t_stat:.4f}")
print(f"p-value:     {p_val:.4e}")

if p_val < 0.05:
    print("Result: Statistically Significant (p < 0.05)")
else:
    print("Result: Not Statistically Significant")


# ## Now plot PGR + pgthAB for catechin degraders
# 

# In[32]:


PGR_only_df = df.loc[(df['MAG'].isin(catechin_degraders)) & (df['enzyme'] == 'PGR')].sort_values(by=['MAG'])
PGR_total_per_sample = PGR_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
PGR_total_per_sample.columns = ['Sample', 'total_getmm']
PGR_total_per_sample = PGR_total_per_sample.merge(replicate_frame, on='Sample', how='left')

# Calculate mean and 95% confidence interval for each treatment-day combination
summary = PGR_total_per_sample.groupby(['treatment', 'day'])['total_getmm'].agg([
    ('mean', 'mean'),
    ('std', 'std'),
    ('count', 'count')
]).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Define colors for treatments
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Plot each treatment
for treatment in summary['treatment'].unique():
    treatment_data = summary[summary['treatment'] == treatment]
    
    # Plot replicate points
    replicate_data = PGR_total_per_sample[PGR_total_per_sample['treatment'] == treatment]
    ax.scatter(replicate_data['day'], replicate_data['total_getmm'],
              color=colors.get(treatment, 'gray'), alpha=0.6, s=50, zorder=3)
    
    # Plot mean line (no markers)
    ax.plot(treatment_data['day'], treatment_data['mean'], 
            linewidth=2.5,
            label=treatment, color=colors.get(treatment, 'gray'))
    
# Formatting
ax.set_xlabel('Day', fontsize=12, fontweight='bold')
ax.set_ylabel('Total getmm', fontsize=12, fontweight='bold')
ax.set_title('PGR Total Over Time by Treatment', 
             fontsize=14, fontweight='bold')
ax.legend(title='Treatment', fontsize=10)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xticks([0, 7, 14, 21, 35])

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10A-PGR_total_over_time_catechin-degraders.pdf', dpi=300)
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(summary.to_string(index=False))


# In[33]:


pgth_only_df = df.loc[(df['MAG'].isin(catechin_degraders)) & (df['enzyme'] == 'pgthAB')].sort_values(by=['MAG'])
pgth_total_per_sample = pgth_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
pgth_total_per_sample.columns = ['Sample', 'total_getmm']
pgth_total_per_sample = pgth_total_per_sample.merge(replicate_frame, on='Sample', how='left')


# Calculate mean and 95% confidence interval for each treatment-day combination
summary = pgth_total_per_sample.groupby(['treatment', 'day'])['total_getmm'].agg([
    ('mean', 'mean'),
    ('std', 'std'),
    ('count', 'count')
]).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Define colors for treatments
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Plot each treatment
for treatment in summary['treatment'].unique():
    treatment_data = summary[summary['treatment'] == treatment]
    
    # Plot replicate points
    replicate_data = pgth_total_per_sample[pgth_total_per_sample['treatment'] == treatment]
    ax.scatter(replicate_data['day'], replicate_data['total_getmm'],
              color=colors.get(treatment, 'gray'), alpha=0.6, s=50, zorder=3)
    
    # Plot mean line (no markers)
    ax.plot(treatment_data['day'], treatment_data['mean'], 
            linewidth=2.5,
            label=treatment, color=colors.get(treatment, 'gray'))
    
# Formatting
ax.set_xlabel('Day', fontsize=12, fontweight='bold')
ax.set_ylabel('Total getmm', fontsize=12, fontweight='bold')
ax.set_title('pgth Total Over Time by Treatment', 
             fontsize=14, fontweight='bold')
ax.legend(title='Treatment', fontsize=10)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xticks([0, 7, 14, 21, 35])

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10B-pgth_total_over_time_catechin-degraders.pdf', dpi=300)
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(summary.to_string(index=False))


# In[34]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# 1. Configuration
# Mapping enzymes and groups to create the 2x2 grid
# Top row: Catechin Degraders | Bottom row: Additional Polyphenol Degraders
plot_configs = [
    {'enzyme': 'PGR',    'mags': catechin_degraders,    'title': 'PGR (Catechin degraders)'},
    {'enzyme': 'pgthAB', 'mags': catechin_degraders,    'title': 'pgthAB (Catechin degraders)'},
    {'enzyme': 'PGR',    'mags': pgthab_mags,           'title': 'PGR (Additional polpyhenol degraders)'},
    {'enzyme': 'pgthAB', 'mags': pgthab_mags,           'title': 'pgthAB (Additional polyphenol degraders)'}
]

colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

def get_sig_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

# 2. Setup Figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()

for i, config in enumerate(plot_configs):
    ax = axes[i]
    
    # Filter data for specific enzyme and specific MAG list
    filtered_df = df.loc[
        (df['enzyme'] == config['enzyme']) & 
        (df['MAG'].isin(config['mags']))
    ].copy()
    
    if filtered_df.empty:
        ax.set_title(f"{config['title']} (No Data)")
        continue

    # Aggregate by sample
    sample_sums = filtered_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
    sample_sums.columns = ['Sample', 'total_getmm']
    sample_sums = sample_sums.merge(replicate_frame, on='Sample', how='left')
    
    # Calculate Welch's T-test (Catechin vs Unamended)
    group1 = sample_sums[sample_sums['treatment'] == 'catechin']['total_getmm']
    group2 = sample_sums[sample_sums['treatment'] == 'unamended']['total_getmm']
    
    sig_text = ""
    if len(group1) > 1 and len(group2) > 1:
        _, p_val = stats.ttest_ind(group1, group2, equal_var=False)
        sig_text = f" ({get_sig_label(p_val)})"

    # Calculate means for line plotting
    summary = sample_sums.groupby(['treatment', 'day'])['total_getmm'].mean().reset_index()

    # 3. Plotting
    for treatment in ['catechin', 'unamended']:
        color = colors.get(treatment, 'gray')
        
        # Replicates (scatter)
        reps = sample_sums[sample_sums['treatment'] == treatment]
        ax.scatter(reps['day'], reps['total_getmm'], color=color, alpha=0.5, s=40, zorder=3)
        
        # Mean Trend (line)
        trend = summary[summary['treatment'] == treatment]
        ax.plot(trend['day'], trend['total_getmm'], color=color, linewidth=2.5, label=treatment, zorder=4)

    # 4. Formatting
    ax.set_title(f"{config['title']}{sig_text}", fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xticks([0, 7, 14, 21, 35])
    
    # Labels (only on edges)
    if i >= 2: ax.set_xlabel('Day', fontweight='bold')
    if i % 2 == 0: ax.set_ylabel('Total geTMM', fontweight='bold')
    
    if i == 0:
        ax.legend(title='Treatment', fontsize=9)

plt.tight_layout()

# 5. Save and Show
save_path = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10_Combined_Enzyme_Panels.pdf'
plt.savefig(save_path, dpi=300)
plt.show()

print(f"Plot saved to: {save_path}")


# In[35]:


# --- 1. Use Your Exact Configuration Lists ---
plot_configs = [
    {'enzyme': 'PGR',    'mags': catechin_degraders,    'functional_group': 'Catechin degraders'},
    {'enzyme': 'pgthAB', 'mags': catechin_degraders,    'functional_group': 'Catechin degraders'},
    {'enzyme': 'PGR',    'mags': pgthab_mags,           'functional_group': 'Additional polyphenol degraders'},
    {'enzyme': 'pgthAB', 'mags': pgthab_mags,           'functional_group': 'Additional polyphenol degraders'}
]

compiled_records = []

# --- 2. Iterate through Panels & Generate Metrics ---
for config in plot_configs:
    # Filter base data for specific enzyme and targeted MAG list
    filtered_df = df.loc[
        (df['enzyme'] == config['enzyme']) & 
        (df['MAG'].isin(config['mags']))
    ].copy()
    
    if filtered_df.empty:
        continue

    # Aggregate total expression by sample
    sample_sums = filtered_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
    sample_sums.columns = ['Sample', 'total_getmm']
    sample_sums = sample_sums.merge(replicate_frame, on='Sample', how='left')
    
    # Calculate Welch's T-test p-value across treatments (collapsing days)
    group1 = sample_sums[sample_sums['treatment'] == 'catechin']['total_getmm']
    group2 = sample_sums[sample_sums['treatment'] == 'unamended']['total_getmm']
    
    p_val_numeric = np.nan
    if len(group1) > 1 and len(group2) > 1:
        _, p_val = stats.ttest_ind(group1, group2, equal_var=False)
        p_val_numeric = round(p_val, 6)

    # Append individual biological replicate rows to compilation array
    for _, row in sample_sums.iterrows():
        compiled_records.append({
            'Panel_Title': f"{config['enzyme']} ({config['functional_group']})",
            'Enzyme': config['enzyme'],
            'Target_Functional_Group': config['functional_group'],
            'Sample_ID': row['Sample'],
            'Treatment': row['treatment'].capitalize(),
            'Timepoint_Days': int(row['day']),
            'Total_Expression_geTMM': round(row['total_getmm'], 4),
            'Welch_Anova_p_value': p_val_numeric
        })

# --- 3. Build tidy DataFrame ---
df_supplemental_s11 = pd.DataFrame(compiled_records)

# Sort logically for spreadsheet readability
df_supplemental_s11 = df_supplemental_s11.sort_values(
    by=['Target_Functional_Group', 'Enzyme', 'Treatment', 'Timepoint_Days', 'Sample_ID']
).reset_index(drop=True)

# --- 4. Export to Destination Directory ---
export_path = '/users/PAS1573/riddell26/data/S11_enzyme_grid_data.csv'
df_supplemental_s11.to_csv(export_path, index=False)
df_supplemental_s11


# In[36]:


MAG_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_GTDB_getmm.tsv', sep='\t')
MAG_metaT.head()


# In[37]:


# Define colors for treatments
treatment_colors = {
    'catechin': '#FC9B2D',
    'unamended': '#7ACBC3'
}

# Load metadata once (outside functions)
metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv',
                  usecols=['Sample', 'treatment', 'timepoint'])
metadata['day'] = metadata['timepoint'].apply(lambda x: int(x.replace('day', '')))
metadata = metadata.drop(columns=['timepoint'])



# In[38]:


df['geTMM'].min()


# In[39]:


# Create a mapping of MAG to taxonomic rank
mag_to_tax = df.groupby('MAG')['highest_host_tax_rank'].first().to_dict()

# Identify MAGs enriched in catechin at different numbers of timepoints (days 14-35)
enriched_3plus_timepoints = set()
enriched_2_timepoints = set()
enriched_1_timepoint = set()
not_enriched = set()

for mag in all_mags:
    enriched_timepoints = 0
    
    # Check both PGR and pgthAB
    for enzyme in ['PGR', 'pgthAB']:
        mag_enzyme_data = df[(df['MAG'] == mag) & (df['enzyme'] == enzyme)]
        
        # Check days 14, 21, 35
        for time in [14, 21, 35]:
            catechin_vals = mag_enzyme_data[(mag_enzyme_data['treatment'] == 'catechin') & 
                                           (mag_enzyme_data['time'] == time)]['geTMM'].values
            unamended_vals = mag_enzyme_data[(mag_enzyme_data['treatment'] == 'unamended') & 
                                            (mag_enzyme_data['time'] == time)]['geTMM'].values
            
            # Check if we have data for both treatments
            if len(catechin_vals) > 0 and len(unamended_vals) > 0:
                max_catechin = catechin_vals.max()
                mean_catechin = catechin_vals.mean()
                max_unamended = unamended_vals.max()
                mean_unamended = unamended_vals.mean()
                
                # Check for enrichment: BOTH max catechin > max unamended AND mean catechin > mean unamended
                if max_catechin > max_unamended and mean_catechin > mean_unamended:
                    enriched_timepoints += 1
    
    # Categorize based on number of enriched timepoints
    if enriched_timepoints >= 3:
        enriched_3plus_timepoints.add(mag)
    elif enriched_timepoints == 2:
        enriched_2_timepoints.add(mag)
    elif enriched_timepoints == 1:
        enriched_1_timepoint.add(mag)
    else:
        not_enriched.add(mag)

# Sort: 3+ timepoints first, then 2, then 1, then not enriched
sorted_mags = (sorted(enriched_3plus_timepoints) + 
               sorted(enriched_2_timepoints) + 
               sorted(enriched_1_timepoint) + 
               sorted(not_enriched))

print(f"Enriched in 3+ timepoints: {len(enriched_3plus_timepoints)}")
print(f"Enriched in 2 timepoints: {len(enriched_2_timepoints)}")
print(f"Enriched in 1 timepoint: {len(enriched_1_timepoint)}")
print(f"Not enriched: {len(not_enriched)}")
print(f"Total: {len(sorted_mags)}")

# Extract taxonomic ranks for each group and count unique values
print("\n=== Taxonomic Rank Distribution ===")

for group_name, group_mags in [
    ("Enriched in 3+ timepoints", enriched_3plus_timepoints),
    ("Enriched in 2 timepoints", enriched_2_timepoints),
    ("Enriched in 1 timepoint", enriched_1_timepoint),
    ("Not enriched", not_enriched)
]:
    # Get taxonomic ranks for MAGs in this group
    tax_ranks = [mag_to_tax.get(mag, 'Unknown') for mag in group_mags]
    unique_ranks = set(tax_ranks)
    
    print(f"\n{group_name}:")
    print(f"  Total MAGs: {len(group_mags)}")
    print(f"  Unique taxonomic ranks: {len(unique_ranks)}")
    print(f"  Ranks: {sorted(unique_ranks)}")


# # Now that we found the enriched MAGs, add the hydrogenase activity

# In[40]:


hydrogenase = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/hydrogenase_expression.csv')
hydrogenase


# In[41]:


# group by uptake, bidirectional, and consuming
hyd_direction = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/hydrogenase_key.csv')
hyd_direction


# In[42]:


hyd_dict = hyd_direction.groupby('direction')['hydrogenase_group'].apply(list).to_dict()
hyd_dict


# In[43]:


hyd_direction.direction.unique()


# In[44]:


hydrogenase = hydrogenase.merge(hyd_direction, on='hydrogenase_group', how='left')


# In[45]:


hydrogenase


# In[46]:


enriched_mags = (sorted(enriched_3plus_timepoints) + 
               sorted(enriched_2_timepoints) + 
               sorted(enriched_1_timepoint)
                )
hydrogenase_enriched_only = hydrogenase.loc[hydrogenase['MAG'].isin(enriched_mags)]
len(enriched_mags), hydrogenase_enriched_only.MAG.nunique()


# In[47]:


hydrogenase_enriched_only = hydrogenase_enriched_only.merge(df[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination']].drop_duplicates(), how='left')


# In[48]:


hydrogenase_enriched_only = hydrogenase_enriched_only.rename(columns={'hydrogenase_group': 'description'})


# In[49]:


h_melted = hydrogenase_enriched_only.melt(id_vars=['gene', 'description', 'direction', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination'],
                              var_name='sample', value_name='geTMM')
h_melted


# In[50]:


# no longer grouping by enzyme group, instead by direction.
# h_melted['enzyme'] = h_melted['description'].apply(lambda x: x.split('_')[0])

h_melted = h_melted.rename(columns={'direction': 'enzyme'})


# In[51]:


h_melted = h_melted.groupby(['gene', 'enzyme', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'sample', 'Completeness', 'Contamination']).agg({"geTMM": "sum"}).reset_index()


# In[52]:


h_melted.columns


# In[53]:


h_melted = h_melted.merge(df[['sample', 'treatment', 'time']].drop_duplicates(), on='sample', how='left')
h_melted.columns


# In[54]:


df = df[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[55]:


h_melted = h_melted[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[56]:


h_melted.MAG.nunique()


# In[57]:


df_merged = pd.concat([h_melted, df])


# In[58]:


df_merged


# In[59]:


# Convert to dictionary with MAG as key and (completeness, contamination) tuple as value
completeness_contamination_dict = dict(zip(
    completeness['MAG'],
    zip(completeness['Completeness'], completeness['Contamination'])
))


# In[60]:


df_merged = df_merged.loc[(df_merged['Completeness'] > 50) & (df_merged['Contamination'] < 10)]
df_merged


# # Add additional polpyhenol metabolisms

# In[61]:


df_merged.columns


# In[62]:


camper_genes_df = camper_genes_df.merge(df[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination']].drop_duplicates(), how='left', on='MAG')
camper_genes_df = camper_genes_df.loc[camper_genes_df['MAG'].isin(enriched_mags)]
camper_genes_df


# In[63]:


camper_genes_df.columns


# In[64]:


camper_genes_melted = camper_genes_df.melt(id_vars=['gene', 'fasta', 'scaffold',
       'gene_position', 'camper_hits', 'camper_rank', 'camper_bitScore',
       'camper_id', 'camper_definition', 'camper_search_type', 'MAG', 'GTDB',
       'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank',
       'Completeness', 'Contamination'], var_name='sample', value_name='geTMM')

camper_genes_melted['enzyme'] = camper_genes_melted['camper_definition'].fillna('').apply(lambda x: x.split(';', 1)[0])


# In[65]:


camper_genes_melted


# In[66]:


camper_genes_melted.enzyme.unique()


# In[67]:


camper_genes_melted.camper_definition.unique()


# In[68]:


# filter for pgthAB and PGR-containing MAGs
# pgthab_pgr_mags = pgr_mags | pgthab_mags
# pgthab_pgr_mags_genes_melted = camper_genes_melted.loc[camper_genes_melted['fasta'].isin(pgthab_pgr_mags)]


# In[69]:


# what enzymes are upstream of phloroglucinol? Find these and put them in a list, then filter.
additional_enzymes = [
 'FLR', 'CHI', 'FCR', 'PHY', # Flavonoids
 'CalB', # Phenyl-propenes
# 'CLD', # Lignans
 'EDL', # Lignans
 'PhcD', # Lignans
 'BER', # Lignans
 'CurA', # Curcuminoids
 'GRAD', # Phenolic acids
 'sam5', # Isoflavonoids, cinnamic acids
 'PhcC', # Lignans
 'ChlE', # Cinnamic acid metabolism
 'FdeE', # Flavonones
 'FdeC', # Flavonones
 'IEM', # Phenyl-propenes
 'SesA', # Lignans
 'PinZ', # Lignans
 'PAH', # Lignans
 'GLM', # Lignans
 'EhyB', # Phenyl-propenes
 'CalA', # Phenyl-propenes
 'FCS', # Non-specific
 'PhcG' # Lignans
 # 'CarA', # Phenolic acids
 # 'CarB', # Phenolic acids
 # 'CarC', # Phenolic acids
 # 'CarD', # Phenolic acids
 # 'CarE' # Phenolic acids
]


# In[70]:


additional_melted = camper_genes_melted.loc[camper_genes_melted['enzyme'].isin(additional_enzymes)]
additional_melted


# In[71]:


additional_melted = additional_melted.merge(df_merged[['sample', 'treatment', 'time']].drop_duplicates(), on='sample', how='left')


# In[72]:


additional_melted.MAG.nunique()


# In[73]:


df_merged.columns


# In[74]:


additional_melted_to_merge = additional_melted[['gene', 'MAG', 'enzyme', 'sample', 'treatment', 'time', 'geTMM', 'GTDB',
       'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank',
       'Completeness', 'Contamination']]


# In[75]:


df_merged.MAG.nunique()


# In[76]:


df_merged = pd.concat([df_merged, additional_melted_to_merge])


# In[77]:


df_merged.enzyme.unique()


# # Make heatmap with df_merged

# In[78]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Define enzymes to include
enzymes = ['PGR', 'pgthAB', 'carABCDE', 'CLD', 'mhpABCD', 'hpaAB', 'hpaDEFGH', 'uptake', 'bifurcating', 'sensing', 'bidirectional', 'producing']

# Define timepoints to include
timepoints = [14, 21, 35]

# Pseudocount to avoid divide-by-zero
pseudocount = 0.0001

def calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus, enriched_2, enriched_1, mag_to_tax, df_merged):
    """
    Calculate log2 fold-change (catechin/unamended) for each MAG, enzyme, and timepoint.
    Only includes MAGs enriched in at least one timepoint.
    
    Returns:
    --------
    fc_df : DataFrame
        Dataframe with columns for each enzyme-timepoint combination
    """
    
    # Combine all enriched MAGs and sort
    sorted_mags = sorted(enriched_3plus) + sorted(enriched_2) + sorted(enriched_1)
    
    # Create column names for each enzyme-timepoint combination
    column_names = []
    for enzyme in enzymes:
        for time in timepoints:
            column_names.append(f"{enzyme}_day{time}")
    
    # Initialize matrix to store fold-changes
    fc_matrix = np.full((len(sorted_mags), len(column_names)), np.nan)
    
    # Calculate fold-change for each MAG, enzyme, and timepoint
    for i, mag in enumerate(sorted_mags):
        col_idx = 0
        for enzyme in enzymes:
            for time in timepoints:
                # Get data for this MAG, enzyme, and timepoint
                mag_enzyme_time_data = df_merged[(df_merged['MAG'] == mag) & 
                                              (df_merged['enzyme'] == enzyme) & 
                                              (df_merged['time'] == time)]
                
                if len(mag_enzyme_time_data) == 0:
                    # No data - leave as NaN (will be black in heatmap)
                    col_idx += 1
                    continue
                
                # Calculate mean expression for each treatment at this timepoint
                catechin_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'catechin']['geTMM'].mean()
                unamended_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'unamended']['geTMM'].mean()
                
                # Handle NaN values
                if pd.isna(catechin_mean):
                    catechin_mean = 0
                if pd.isna(unamended_mean):
                    unamended_mean = 0
                
                # Calculate log2 fold-change with pseudocount
                log2_fc = np.log2((catechin_mean + pseudocount) / (unamended_mean + pseudocount))
                fc_matrix[i, col_idx] = log2_fc
                
                col_idx += 1
    
    # Create DataFrame
    fc_df = pd.DataFrame(fc_matrix, index=sorted_mags, columns=column_names)
    
    # Add taxonomic information
    fc_df['Taxonomy'] = sorted_mags
    
    # Add enrichment type
    enrichment_type = []
    for mag in sorted_mags:
        if mag in enriched_3plus:
            enrichment_type.append('3 timepoints')
        elif mag in enriched_2:
            enrichment_type.append('2 timepoints')
        elif mag in enriched_1:
            enrichment_type.append('1 timepoint')
    fc_df['Enrichment_Type'] = enrichment_type
    
    return fc_df

def plot_foldchange_heatmap_enriched_only(fc_df, mag_to_tax):
    """
    Create heatmap of log2 fold-changes with timepoints for each enzyme.
    Enzymes on y-axis, MAGs on x-axis, with timepoints as sub-columns.
    Only shows enriched MAGs.
    """
    
    # Separate annotation columns from data
    data_cols = [col for col in fc_df.columns if col not in ['Taxonomy', 'Enrichment_Type']]
    heatmap_data = fc_df[data_cols]
    
    # Transpose the data (enzymes will be rows)
    heatmap_data = heatmap_data.T
    
    # Split data by enrichment type for column organization
    enriched_3plus_idx = fc_df[fc_df['Enrichment_Type'] == '3 timepoints'].index
    enriched_2_idx = fc_df[fc_df['Enrichment_Type'] == '2 timepoints'].index
    enriched_1_idx = fc_df[fc_df['Enrichment_Type'] == '1 timepoint'].index
    
    # Reorder columns by enrichment type
    column_order = list(enriched_3plus_idx) + list(enriched_2_idx) + list(enriched_1_idx)
    heatmap_data = heatmap_data[column_order]
    
    # Get taxonomy labels and enrichment types in same order
    taxonomy_labels = [f"{mag}\n{mag_to_tax.get(mag, 'Unknown')}; {completeness_contamination_dict.get(mag, 'Unknown')}" for mag in column_order]
    enrichment_types = fc_df.loc[column_order, 'Enrichment_Type'].values
    
    # Create figure
    fig = plt.figure(figsize=(max(16, len(column_order) * 0.5), 10))
    ax = fig.add_subplot(111)
    
    # Define color map - use a colormap that doesn't include white
    # RdBu_r goes from red (positive) through light colors to blue (negative)
    # We'll use a different palette that avoids white in the middle
    from matplotlib.colors import LinearSegmentedColormap
    
    # Create custom colormap: dark blue -> light blue -> light red -> dark red
    # This avoids pure white in the data range
    colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
              '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
    cmap.set_bad(color='#e0e0e0')  # Light gray for missing data
    
    # Determine vmin and vmax for symmetric color scale
    vmax = np.nanmax(np.abs(heatmap_data.values))
    vmin = -vmax
    
    # Create custom row labels (just enzyme names, timepoints will be in the cells)
    row_labels = []
    for col in heatmap_data.index:
        enzyme = col.split('_')[0]
        day = col.split('_day')[1]
        row_labels.append(f"{enzyme}_{day}")
    
    # Plot heatmap
    sns.heatmap(heatmap_data, 
                cmap=cmap,
                center=0,
                vmin=vmin,
                vmax=vmax,
                ax=ax,
                cbar=False,
                yticklabels=row_labels,
                xticklabels=taxonomy_labels,
                linewidths=0.5,
                linecolor='lightgray')
    
    # Add thicker black lines to separate enzymes (every 4 rows = 4 timepoints)
    for i in range(len(enzymes)):
        if i > 0:  # Don't draw line before first enzyme
            ax.axhline(y=i*len(timepoints), color='black', linewidth=3)
    
    # Add thicker vertical lines to separate MAGs
    for i in range(len(column_order)):
        if i > 0:
            ax.axvline(x=i, color='black', linewidth=1.5)
    
    ax.set_ylabel('Enzyme_Day', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9, length=0)
    
    # Set x-axis labels (taxonomy) at top and rotate 45 degrees
    ax.tick_params(axis='x', labelsize=7, rotation=45, labelbottom=False, labeltop=True, 
                  bottom=False, top=True)
    ax.xaxis.set_label_position('top')
    
    # Adjust label alignment for readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    
    # Add white vertical lines to separate enrichment groups
    boundaries = []
    prev_type = None
    for i, enrich_type in enumerate(enrichment_types):
        if enrich_type != prev_type and prev_type is not None:
            boundaries.append(i)
        prev_type = enrich_type
    
    for boundary in boundaries:
        ax.axvline(x=boundary, color='white', linewidth=4, linestyle='-')
    
    # Add text labels for each group at the top
    group_starts = [0] + boundaries
    group_ends = boundaries + [len(column_order)]
    
    for start, end in zip(group_starts, group_ends):
        midpoint = (start + end) / 2
        enrich_type = enrichment_types[start]
        ax.text(midpoint, len(heatmap_data) + 1.2, f'Catechin-enriched\n{enrich_type}',
               va='bottom', ha='center', fontsize=10, fontweight='bold',
               rotation=0)
    
    # Add colorbar at bottom
    cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.015])
    
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('log2-fold-change catechin / unamended mean GeTMM', fontsize=11, fontweight='bold')
    cbar.ax.tick_params(labelsize=9)
    
    plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_enzyme_response_heatmap_enriched_only_hyd_direction.pdf',
                dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig

# Generate the heatmap with only enriched MAGs
fc_df_enriched = calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus_timepoints, 
                                                                  enriched_2_timepoints, 
                                                                  enriched_1_timepoint,
                                                                  mag_to_tax, df_merged)
fig = plot_foldchange_heatmap_enriched_only(fc_df_enriched, mag_to_tax)

# Display summary statistics
print(f"\nSummary Statistics (Enriched MAGs Only):")
print(f"Total MAGs shown: {len(fc_df_enriched)}")
print(f"Enriched 3+ timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '3 timepoints')}")
print(f"Enriched 2 timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '2 timepoints')}")
print(f"Enriched 1 timepoint: {sum(fc_df_enriched['Enrichment_Type'] == '1 timepoint')}")
print(f"\nLog2 FC range: {np.nanmin(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f} to {np.nanmax(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f}")


# In[79]:


# --- 1. Reshape the Heatmap Data into a Tidy Format ---
# (Assuming fc_df_enriched from your heatmap script is already generated)
# Reset index to turn the MAG identifier into a formal column
supplemental_5a = fc_df_enriched.reset_index().rename(columns={'index': 'MAG_ID'})

# Identify data columns (enzyme_day combinations) vs. annotation columns
annotation_cols = ['MAG_ID', 'Taxonomy', 'Enrichment_Type']
data_cols = [col for col in supplemental_5a.columns if col not in annotation_cols]

# Melt the DataFrame from wide format to a long/tidy format
df_tidy_5a = supplemental_5a.melt(
    id_vars=annotation_cols,
    value_vars=data_cols,
    var_name='Enzyme_Timepoint',
    value_name='Log2_Fold_Change'
)

# --- 2. Clean and Parse Merged Fields for Reader Clarity ---
# Split the "Enzyme_Timepoint" column into two explicit scientific columns
df_tidy_5a['Enzyme'] = df_tidy_5a['Enzyme_Timepoint'].apply(lambda x: x.split('_day')[0])
df_tidy_5a['Timepoint_Days'] = df_tidy_5a['Enzyme_Timepoint'].apply(lambda x: int(x.split('_day')[1]))

# Extract clean taxonomic lineages and quality values using your dictionaries
# (Avoids relying on the plot-formatted label strings containing '\n' characters)
df_tidy_5a['Taxonomic_Lineage'] = df_tidy_5a['MAG_ID'].map(lambda m: mag_to_tax.get(m, 'Unknown'))
df_tidy_5a['Completeness_Contamination'] = df_tidy_5a['MAG_ID'].map(lambda m: completeness_contamination_dict.get(m, 'Unknown'))

# Organize column order logically for publication spreadsheet guidelines
final_column_order = [
    'MAG_ID',
    'Taxonomic_Lineage',
    'Completeness_Contamination',
    'Enrichment_Type',
    'Enzyme',
    'Timepoint_Days',
    'Log2_Fold_Change'
]
df_tidy_5a = df_tidy_5a[final_column_order]

# Sort columns systematically so it reads intuitively
df_tidy_5a = df_tidy_5a.sort_values(
    by=['Enrichment_Type', 'MAG_ID', 'Enzyme', 'Timepoint_Days']
).reset_index(drop=True)

# --- 3. Format Precision and Handle Missing Values ---
# Round numeric entries to a publication-standard precision
df_tidy_5a['Log2_Fold_Change'] = df_tidy_5a['Log2_Fold_Change'].round(4)

# Replace internal NaN representations with an explicit descriptor for missing data
# This represents cases where there was no sequence coverage (rendered as light gray boxes in the plot)
df_tidy_5a['Log2_Fold_Change'] = df_tidy_5a['Log2_Fold_Change'].fillna('No Detected Expression')

# --- 4. Export to Journal-Ready CSV ---
csv_path_5a = '/users/PAS1573/riddell26/data/5A_MAG_polyphenol_activity_heatmap_source_data.csv'
df_tidy_5a.to_csv(csv_path_5a, index=False)
df_tidy_5a


# In[80]:


import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats  # Added for t-test

# Define the directions
hydro_directions = ['uptake', 'bifurcating', 'sensing', 'producing']
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

def get_sig_label(p):
    """Helper to convert p-value to significance stars."""
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

# Create the 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()

for i, direction in enumerate(hydro_directions):
    ax = axes[i]
    
    # 1. Filter using df_merged
    dir_only_df = df_merged.loc[(df_merged['enzyme'] == direction)].copy()
    
    if dir_only_df.empty:
        ax.set_title(f'{direction.capitalize()} (No Data)')
        continue

    # 2. Aggregation and Merge
    dir_total_per_sample = dir_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
    dir_total_per_sample.columns = ['Sample', 'total_getmm']
    dir_total_per_sample = dir_total_per_sample.merge(replicate_frame, on='Sample', how='left')

    time_col = 'day' if 'day' in dir_total_per_sample.columns else 'time'
    
    # --- 3. Perform Welch's T-Test ---
    # We compare the 'total_getmm' values between the two treatments
    group1 = dir_total_per_sample[dir_total_per_sample['treatment'] == 'catechin']['total_getmm']
    group2 = dir_total_per_sample[dir_total_per_sample['treatment'] == 'unamended']['total_getmm']
    
    sig_text = ""
    if not group1.empty and not group2.empty:
        t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=False)
        sig_text = f" ({get_sig_label(p_val)})"

    # 4. Calculate mean stats for plotting lines
    summary = dir_total_per_sample.groupby(['treatment', time_col])['total_getmm'].agg([
        ('mean', 'mean')
    ]).reset_index()

    # 5. Plotting
    for treatment in summary['treatment'].unique():
        t_summary = summary[summary['treatment'] == treatment]
        replicates = dir_total_per_sample[dir_total_per_sample['treatment'] == treatment]
        
        color = colors.get(treatment, 'gray')
        
        # Plot replicates
        ax.scatter(replicates[time_col], replicates['total_getmm'],
                   color=color, alpha=0.5, s=30, zorder=3)
        
        # Plot mean line
        ax.plot(t_summary[time_col], t_summary['mean'], 
                linewidth=2, label=treatment, color=color, zorder=4)
        
    # 6. Facet Formatting
    # Title now includes the significance label
    ax.set_title(f'{direction.capitalize()}{sig_text}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xticks([0, 7, 14, 21, 35])
    
    if i >= 2: ax.set_xlabel('Day', fontweight='bold')
    if i % 2 == 0: ax.set_ylabel('Total geTMM', fontweight='bold')
    
    if i == 0:
        ax.legend(title='Treatment', fontsize=9)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10B-Hydrogenase_Facet_Grid.pdf', dpi=300)
plt.show()


# In[81]:


catechin_MAG_df = hydrogenase.loc[hydrogenase['MAG'].isin(catechin_degraders)]
catechin_MAG_df.columns


# In[82]:


import matplotlib.pyplot as plt
import pandas as pd

# 1. Reshape catechin_MAG_df from Wide to Long format
# Identify the sample columns (those starting with 'STM')
sample_cols = [col for col in catechin_MAG_df.columns if col.startswith('STM')]

df_long = catechin_MAG_df.melt(
    id_vars=['gene', 'MAG', 'direction'], 
    value_vars=sample_cols,
    var_name='Sample', 
    value_name='geTMM'
)

# 2. Merge with replicate_frame to get 'treatment' and 'day'
# This creates the 'df_merged' equivalent needed for the plot
df_plot_ready = df_long.merge(replicate_frame[['Sample', 'treatment', 'day']], on='Sample', how='left')

# 3. Plotting logic using the reshaped dataframe
hydro_directions = ['uptake', 'bifurcating', 'sensing', 'producing']
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()

for i, direction in enumerate(hydro_directions):
    ax = axes[i]
    
    # Filter for the specific functional direction
    dir_only_df = df_plot_ready.loc[df_plot_ready['direction'] == direction].copy()
    
    if dir_only_df.empty:
        ax.set_title(f'Hydrogenase: {direction.capitalize()} (No Data)')
        continue

    # Aggregation: Sum geTMM per sample within this hydrogenase direction
    dir_total_per_sample = dir_only_df.groupby(['Sample', 'treatment', 'day']).agg({'geTMM': 'sum'}).reset_index()

    # Calculate mean stats for the lines
    summary = dir_total_per_sample.groupby(['treatment', 'day'])['geTMM'].agg([
        ('mean', 'mean')
    ]).reset_index()

    # Plotting
    for treatment in summary['treatment'].unique():
        t_summary = summary[summary['treatment'] == treatment]
        replicates = dir_total_per_sample[dir_total_per_sample['treatment'] == treatment]
        
        color = colors.get(treatment, 'gray')
        
        # Plot replicate points
        ax.scatter(replicates['day'], replicates['geTMM'],
                   color=color, alpha=0.5, s=30, zorder=3)
        
        # Plot mean line
        ax.plot(t_summary['day'], t_summary['mean'], 
                linewidth=2, label=treatment, color=color, zorder=4)
        
    # Formatting
    ax.set_title(f'Hydrogenase: {direction.capitalize()}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xticks([0, 7, 14, 21, 35])
    
    if i >= 2: ax.set_xlabel('Day', fontweight='bold')
    if i % 2 == 0: ax.set_ylabel('Total geTMM', fontweight='bold')
    if i == 0: ax.legend(title='Treatment', fontsize=9)

plt.tight_layout()
plt.show()


# In[83]:


# Plot hydrogenase uptake over time as proportion by catechin degraders vs non-catechin-degraders


# In[84]:


hydrogenase_uptake_catechin = (
    hydrogenase
    .loc[
        (hydrogenase['MAG'].isin(catechin_degraders))
        & (hydrogenase['direction'] == 'uptake')
    ]
)

hydrogenase_uptake_non_degraders = (
    hydrogenase
    .loc[
        (hydrogenase['MAG'].isin(enriched_mags))
        & (hydrogenase['direction'] == 'uptake')
    ]
)

# 1. Helper function to transform data from wide to long format
def transform_to_plot_ready(df, replicate_df):
    # Identify sample columns (those starting with 'STM')
    sample_cols = [col for col in df.columns if col.startswith('STM')]
    
    # Reshape from wide to long
    df_long = df.melt(
        id_vars=['gene', 'MAG', 'direction'], 
        value_vars=sample_cols,
        var_name='Sample', 
        value_name='geTMM'
    )
    
    # Merge with replicate metadata to get treatment and day info
    return df_long.merge(replicate_df[['Sample', 'treatment', 'day']], on='Sample', how='left')

# 2. Process both dataframes
plot_ready_catechin = transform_to_plot_ready(hydrogenase_uptake_catechin, replicate_frame)
plot_ready_non_degraders = transform_to_plot_ready(hydrogenase_uptake_non_degraders, replicate_frame)

# 3. Plotting Configuration
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}
fig, ax = plt.subplots(figsize=(10, 6))

# Groups to iterate through: (dataframe, linestyle, label_suffix)
groups_to_plot = [
    (plot_ready_catechin, '-', 'Degraders'),
    (plot_ready_non_degraders, '--', 'Non-Degraders')
]

# 4. Iterative Plotting Logic
for df_plot_ready, lstyle, suffix in groups_to_plot:
    # Aggregation: Sum geTMM per sample within this subset
    total_per_sample = df_plot_ready.groupby(['Sample', 'treatment', 'day']).agg({'geTMM': 'sum'}).reset_index()

    # Calculate mean stats for the lines
    summary = total_per_sample.groupby(['treatment', 'day'])['geTMM'].agg([('mean', 'mean')]).reset_index()

    # Plot Treatment Lines
    for treatment in summary['treatment'].unique():
        t_summary = summary[summary['treatment'] == treatment]
        replicates = total_per_sample[total_per_sample['treatment'] == treatment]
        color = colors.get(treatment, 'gray')
        
        # Plot individual replicate points
        ax.scatter(replicates['day'], replicates['geTMM'],
                   color=color, alpha=0.4, s=30, zorder=3)
        
        # Plot mean line with specific linestyle (solid vs dashed)
        ax.plot(t_summary['day'], t_summary['mean'], 
                linewidth=2.5, linestyle=lstyle, color=color, 
                label=f"{treatment.capitalize()} ({suffix})", zorder=4)

# 5. Formatting and Aesthetics
ax.set_title('Uptake Hydrogenase expression', fontsize=13, fontweight='bold')
ax.set_xlabel('Day', fontweight='bold')
ax.set_ylabel('Total geTMM', fontweight='bold')
ax.set_xticks([0, 7, 14, 21, 35])
ax.grid(True, alpha=0.2, linestyle='--')

# Place legend outside the plot area
ax.legend(title='Treatment & Group', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/hydrogenase_uptake_cat_vs_additional.pdf', dpi=300)
plt.show()


# In[85]:


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# 1. Processing and Aggregation (assuming previous steps are done)
plot_ready_catechin['subset'] = 'Degraders'
plot_ready_non_degraders['subset'] = 'Non-Degraders'
df_combined = pd.concat([plot_ready_catechin, plot_ready_non_degraders])

df_summed = (
    df_combined.groupby(['Sample', 'treatment', 'day', 'subset'])
    .agg({'geTMM': 'sum'})
    .reset_index()
)
df_summed['group'] = df_summed['treatment'] + " (" + df_summed['subset'] + ")"
df_summed['log_geTMM'] = np.log1p(df_summed['geTMM'])

# 2. Plotting Configuration
group_colors = {
    'catechin (Degraders)': '#E65100',
    'catechin (Non-Degraders)': '#FFB74D',
    'unamended (Degraders)': '#006064',
    'unamended (Non-Degraders)': '#80CBC4'
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1.5, 1]})

# --- Subplot 1: Line Plot ---
linestyles = {'Degraders': '-', 'Non-Degraders': '--'}
for (subset, treatment), group_df in df_summed.groupby(['subset', 'treatment']):
    label_str = f"{treatment} ({subset})"
    color = group_colors.get(label_str, 'gray')
    lstyle = linestyles.get(subset, '-')
    
    stats_summary = group_df.groupby('day')['geTMM'].agg(['mean', 'sem']).reset_index()
    stats_summary['ci_upper'] = stats_summary['mean'] + (1.96 * stats_summary['sem'])
    stats_summary['ci_lower'] = stats_summary['mean'] - (1.96 * stats_summary['sem'])

    ax1.plot(stats_summary['day'], stats_summary['mean'], 
             label=label_str.capitalize(),
             color=color, linestyle=lstyle, linewidth=2.5, zorder=4)
    ax1.fill_between(stats_summary['day'], stats_summary['ci_lower'], stats_summary['ci_upper'], 
                     color=color, alpha=0.15, zorder=3)
    ax1.scatter(group_df['day'], group_df['geTMM'], color=color, alpha=0.3, s=25, zorder=2)

ax1.set_title('Hydrogenase Expression Over Time (Linear Scale)', fontsize=14)
ax1.set_ylabel('Total geTMM')
ax1.set_xticks([0, 7, 14, 21, 35])
ax1.legend(loc='upper left')

# --- Subplot 2: Boxplots (Fixed Warnings) ---
order = ['catechin (Degraders)', 'unamended (Degraders)', 'catechin (Non-Degraders)', 'unamended (Non-Degraders)']

# Added hue=group and legend=False to resolve palette FutureWarning
sns.boxplot(data=df_summed, x='group', y='log_geTMM', order=order, ax=ax2, 
            hue='group', palette=group_colors, legend=False, fliersize=0)
sns.stripplot(data=df_summed, x='group', y='log_geTMM', order=order, ax=ax2, 
              color='black', alpha=0.4, jitter=True)

# Using raw strings (r'') to resolve SyntaxWarning for LaTeX backslashes
ax2.set_title(r'Expression Distribution ($\log(x+1)$ Scale)', fontsize=14)
ax2.set_ylabel(r'$\log(\text{Total geTMM} + 1)$')
ax2.set_xlabel('')

# Fixed UserWarning: Use set_ticks before set_ticklabels
ax2.set_xticks(range(len(order)))
ax2.set_xticklabels([label.replace(' ', '\n') for label in order])

# --- 3. Welch's T-Test ---
cat_deg_log = df_summed[df_summed['group'] == 'catechin (Degraders)']['log_geTMM']
unm_deg_log = df_summed[df_summed['group'] == 'catechin (Non-Degraders)']['log_geTMM']
t_stat, p_val = stats.ttest_ind(cat_deg_log, unm_deg_log, equal_var=False)

y_max = df_summed['log_geTMM'].max()
ax2.text(0.5, y_max * 0.9, f"Welch's p = {p_val:.4f}", 
         ha='center', va='bottom', fontsize=10, color='black')

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/hydrogenase_uptake_cat_vs_additional.pdf', dpi=300)
plt.show()


# In[86]:


# 1. Clean and organize the dataframe for publication
# (Assuming df_summed from your plotting script is already generated)
supplemental_df = df_summed.copy()

# Rename columns to standard, reader-friendly scientific nomenclature
rename_dict = {
    'Sample': 'Sample_ID',
    'treatment': 'Treatment',
    'day': 'Timepoint_Days',
    'subset': 'Functional_Group',
    'geTMM': 'Total_geTMM',
    'log_geTMM': 'Log_Transformed_geTMM_ln_x_plus_1'
}
supplemental_df = supplemental_df.rename(columns=rename_dict)

# Drop columns created purely for matplotlib/seaborn configurations
columns_to_keep = list(rename_dict.values())
supplemental_df = supplemental_df[columns_to_keep]

# Sort logically by group, treatment, and timeline so reviewers/readers can follow it easily
supplemental_df = supplemental_df.sort_values(
    by=['Functional_Group', 'Treatment', 'Timepoint_Days', 'Sample_ID']
).reset_index(drop=True)

# 2. Round floating-point columns to a reasonable precision (e.g., 4 decimal places)
supplemental_df['Total_geTMM'] = supplemental_df['Total_geTMM'].round(4)
supplemental_df['Log_Transformed_geTMM_ln_x_plus_1'] = supplemental_df['Log_Transformed_geTMM_ln_x_plus_1'].round(4)

# 3. Export to Journal-Ready Formats
# Option A: CSV format (Highly preferred by open-science journals)
csv_path = '/users/PAS1573/riddell26/data/5B_hydrogenase_geTMM.csv'
supplemental_df.to_csv(csv_path, index=False)

# Preview the clean table layout
supplemental_df.head()


# In[87]:


# average intensity per mag:

import matplotlib.pyplot as plt
import pandas as pd

# 1. Helper function to transform data (same as previous)
def transform_to_plot_ready(df, replicate_df):
    sample_cols = [col for col in df.columns if col.startswith('STM')]
    df_long = df.melt(
        id_vars=['gene', 'MAG', 'direction'], 
        value_vars=sample_cols,
        var_name='Sample', 
        value_name='geTMM'
    )
    return df_long.merge(replicate_df[['Sample', 'treatment', 'day']], on='Sample', how='left')

# 2. Process dataframes
plot_ready_catechin = transform_to_plot_ready(hydrogenase_uptake_catechin, replicate_frame)
plot_ready_non_degraders = transform_to_plot_ready(hydrogenase_uptake_non_degraders, replicate_frame)

# 3. Plotting Configuration
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}
fig, ax = plt.subplots(figsize=(10, 6))

groups_to_plot = [
    (plot_ready_catechin, '-', 'Degraders'),
    (plot_ready_non_degraders, '--', 'Non-Degraders')
]

# 4. Iterative Plotting Logic (Normalized by MAG count)
for df_plot_ready, lstyle, suffix in groups_to_plot:
    
    # --- CHANGE: Calculate Average per MAG ---
    # First, get the total geTMM and the count of unique MAGs per Sample/Treatment/Day
    sample_stats = df_plot_ready.groupby(['Sample', 'treatment', 'day']).agg(
        total_geTMM=('geTMM', 'sum'),
        mag_count=('MAG', 'nunique')
    ).reset_index()
    
    # Calculate the average geTMM per MAG for each sample
    sample_stats['avg_per_mag'] = sample_stats['total_geTMM'] / sample_stats['mag_count']

    # Calculate mean of these averages for the trend lines
    summary = sample_stats.groupby(['treatment', 'day'])['avg_per_mag'].agg([('mean', 'mean')]).reset_index()

    # Plot Treatment Lines
    for treatment in summary['treatment'].unique():
        t_summary = summary[summary['treatment'] == treatment]
        replicates = sample_stats[sample_stats['treatment'] == treatment]
        color = colors.get(treatment, 'gray')
        
        # Plot individual replicate points (average per MAG)
        ax.scatter(replicates['day'], replicates['avg_per_mag'],
                   color=color, alpha=0.4, s=30, zorder=3)
        
        # Plot mean line
        ax.plot(t_summary['day'], t_summary['mean'], 
                linewidth=2.5, linestyle=lstyle, color=color, 
                label=f"{treatment.capitalize()} ({suffix})", zorder=4)

# 5. Formatting
ax.set_title('Uptake Hydrogenase: Average geTMM per MAG', fontsize=13, fontweight='bold')
ax.set_xlabel('Day', fontweight='bold')
ax.set_ylabel('Average geTMM / MAG', fontweight='bold') # Updated label
ax.set_xticks([0, 7, 14, 21, 35])
ax.grid(True, alpha=0.2, linestyle='--')

ax.legend(title='Treatment & Group', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()


# In[88]:


import matplotlib.pyplot as plt
import pandas as pd

# --- 1. Calculate Totals for Catechin Degraders ---
total_catechin_degraders = plot_ready_catechin['geTMM'].sum()

# --- 2. Split Non-Degraders into Specific Genera vs. Others ---

jaxai_mask = plot_ready_non_degraders['MAG'].str.contains('STM_0716_E_M_E026_E030_E034_A_bin.12', na=False)
jaatfl_mask = plot_ready_non_degraders['MAG'].str.contains('STM_0716_E_M_E058_A_bin.18', na=False)

# Get sums for the specific non-degraders
jaxai_val = plot_ready_non_degraders.loc[jaxai_mask, 'geTMM'].sum()
jaatfl_val = plot_ready_non_degraders.loc[jaatfl_mask, 'geTMM'].sum()

# "Other" Non-Degraders (Total non-degraders minus the two specific ones)
other_non_degraders_val = plot_ready_non_degraders.loc[~(jaxai_mask | jaatfl_mask), 'geTMM'].sum()

# --- 3. Prepare Data for Pie Chart ---
labels = [
    'Catechin degraders', 
    'Actinomycetes g__JAEXAI01\n(polyphenol degrader)', 
    'Actinomycetes g__JAATFL01\n(polyphenol degrader)', 
    'Other polyphenol degraders'
]
sizes = [total_catechin_degraders, jaxai_val, jaatfl_val, other_non_degraders_val]

# Filter out categories with 0 values
plot_data = [(l, s) for l, s in zip(labels, sizes) if s > 0]
labels, sizes = zip(*plot_data)

# --- 4. Plotting ---
# Colors: Shades of grey
colors_pie = ['#FFE0B2', '#FFB74D', '#FF9800', '#E65100']

fig, ax = plt.subplots(figsize=(10, 8))

# Create the pie chart
wedges, texts, autotexts = ax.pie(
    sizes, 
    labels=labels, 
    autopct='%1.1f%%', 
    startangle=140, 
    colors=colors_pie,
    pctdistance=0.75,
    explode=[0.05] * len(labels),
    textprops={'fontsize': 16, 'fontweight': 'bold'}
)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/04-B_hydrogenase_uptake.pdf', dpi=300)
plt.show()


# In[89]:


# --- 1. Calculate Totals for Catechin Degraders ---
total_catechin_degraders = plot_ready_catechin['geTMM'].sum()

# --- 2. Split Non-Degraders into Specific Genera vs. Others ---
jaxai_mask = plot_ready_non_degraders['MAG'].str.contains('STM_0716_E_M_E026_E030_E034_A_bin.12', na=False)
jaatfl_mask = plot_ready_non_degraders['MAG'].str.contains('STM_0716_E_M_E058_A_bin.18', na=False)

jaxai_val = plot_ready_non_degraders.loc[jaxai_mask, 'geTMM'].sum()
jaatfl_val = plot_ready_non_degraders.loc[jaatfl_mask, 'geTMM'].sum()
other_non_degraders_val = plot_ready_non_degraders.loc[~(jaxai_mask | jaatfl_mask), 'geTMM'].sum()

# --- 3. Prepare Tidy Data for Supplemental Table ---
# Reformat labels to be clean and single-lined for spreadsheet cells, 
# and explicitly document the respective target MAG bins for reproducibility.
data_records = [
    {
        'Pie_Chart_Label': 'Catechin degraders',
        'Functional_Group': 'Catechin degraders',
        'Total_Expression_geTMM': total_catechin_degraders
    },
    {
        'Pie_Chart_Label': 'Actinomycetes g__JAEXAI01 (polyphenol degrader)',
        'Functional_Group': 'Polyphenol degraders',
        'Total_Expression_geTMM': jaxai_val
    },
    {
        'Pie_Chart_Label': 'Actinomycetes g__JAATFL01 (polyphenol degrader)',
        'Functional_Group': 'Polyphenol degraders',
        'Total_Expression_geTMM': jaatfl_val
    },
    {
        'Pie_Chart_Label': 'Other polyphenol degraders',
        'Functional_Group': 'Polyphenol degraders',
        'Total_Expression_geTMM': other_non_degraders_val
    }
]

# Convert records to a DataFrame
df_supplemental_pie = pd.DataFrame(data_records)

# Filter out categories with 0 values, matching the exact filtering logic used for the pie chart
df_supplemental_pie = df_supplemental_pie[df_supplemental_pie['Total_Expression_geTMM'] > 0].reset_index(drop=True)

# Pre-calculate percentages so readers can easily map rows back to the pie slices
total_expression = df_supplemental_pie['Total_Expression_geTMM'].sum()
df_supplemental_pie['Slice_Percentage'] = (df_supplemental_pie['Total_Expression_geTMM'] / total_expression) * 100

# Format numerical precision for professional standard publication tables
df_supplemental_pie['Total_Expression_geTMM'] = df_supplemental_pie['Total_Expression_geTMM'].round(4)
df_supplemental_pie['Slice_Percentage'] = df_supplemental_pie['Slice_Percentage'].round(2)

# --- 4. Export to Journal-Ready CSV ---
csv_path = '/users/PAS1573/riddell26/data/5B_hydrogenase_piechart.csv'
df_supplemental_pie.to_csv(csv_path, index=False)

print(f"Supplemental data table successfully saved to: {csv_path}")
print("\nPreview of the supplemental table layout:")
df_supplemental_pie


# In[90]:


import pandas as pd

def get_total_abundance(df):
    """Sums all geTMM values across all sample columns."""
    sample_cols = [col for col in df.columns if col.startswith('STM')]
    return df[sample_cols].sum().sum()

def get_per_gene_abundance(df):
    """Sums geTMM values across samples, grouped by gene."""
    sample_cols = [col for col in df.columns if col.startswith('STM')]
    # Sum across the sample columns for each row, then group by gene name
    gene_sums = df.copy()
    gene_sums['total_abundance'] = df[sample_cols].sum(axis=1)
    return gene_sums.groupby('gene')['total_abundance'].sum()

# 1. Compute Global Abundance and Ratio
total_cat = get_total_abundance(hydrogenase_uptake_catechin)
total_non = get_total_abundance(hydrogenase_uptake_non_degraders)
global_ratio = total_non / total_cat

# 2. Compute Per-Gene Abundance and Ratio
gene_cat = get_per_gene_abundance(hydrogenase_uptake_catechin)
gene_non = get_per_gene_abundance(hydrogenase_uptake_non_degraders)

# Combine into a single comparison dataframe
comparison_df = pd.DataFrame({
    'Catechin_Degraders': gene_cat,
    'Non_Degraders': gene_non
}).fillna(0) # Fill 0 if a gene exists in one group but not the other

comparison_df['Ratio_Non_to_Cat'] = comparison_df['Non_Degraders'] / comparison_df['Catechin_Degraders']

# 3. Display Results
print("--- Global H2 Uptake Abundance Comparison ---")
print(f"Total geTMM (Catechin Degraders): {total_cat:,.2f}")
print(f"Total geTMM (Non-Degraders):      {total_non:,.2f}")
print(f"Non-Degraders are {global_ratio:.2f}x more abundant in total H2 uptake.\n")


# In[91]:


# What proportion of total H2 uptake is attributed to these top two MAGs?
import pandas as pd

def calculate_expression_share(cat_df, non_df, target_labels, replicate_df, tax_df):
    # 1. Combine and transform both datasets to long format
    combined_raw = pd.concat([cat_df, non_df])
    sample_cols = [col for col in combined_raw.columns if col.startswith('STM')]
    
    df_long = combined_raw.melt(
        id_vars=['MAG'], 
        value_vars=sample_cols,
        var_name='Sample', 
        value_name='geTMM'
    )
    
    # 2. Merge taxonomy and create the unique identifier string
    tax_subset = tax_df[['MAG', 'highest_host_tax_rank']].drop_duplicates()
    df_merged = df_long.merge(tax_subset, on='MAG', how='left')
    df_merged['unique_label'] = df_merged['MAG'].astype(str) + " | " + df_merged['highest_host_tax_rank'].astype(str)
    
    # 3. Filter for Catechin treatment and Days 7-35
    df_merged = df_merged.merge(replicate_df[['Sample', 'treatment', 'day']], on='Sample', how='left')
    subset = df_merged[(df_merged['treatment'] == 'catechin') & (df_merged['day'] >= 7)].copy()
    
    # 4. Calculate Total Community Expression
    total_community_uptake = subset['geTMM'].sum()
    
    # 5. Calculate Share for specific target labels
    results = {}
    for label in target_labels:
        label_expression = subset[subset['unique_label'] == label]['geTMM'].sum()
        percentage = (label_expression / total_community_uptake) * 100
        results[label] = {
            'sum_geTMM': label_expression,
            'percentage': percentage
        }
    
    return results, total_community_uptake

# Define the targets
targets = [
    'STM_0716_E_M_E026_E030_E034_A_bin.12 | g__JAEXAI01',
    'STM_0716_E_M_E058_A_bin.18 | g__JAATFL01'
]

# Execute calculation
stats, grand_total = calculate_expression_share(
    hydrogenase_uptake_catechin, 
    hydrogenase_uptake_non_degraders, 
    targets, 
    replicate_frame, 
    MAG_metaT
)

# Output results
print(f"Total H2 Uptake Expression (Days 7-35 Catechin): {grand_total:,.2f} geTMM\n")
for label, data in stats.items():
    print(f"Target: {label}")
    print(f"  - Expression: {data['sum_geTMM']:,.2f} geTMM")
    print(f"  - Attribution: {data['percentage']:.2f}%\n")


# In[92]:


def get_top_mags_plot_data(df, replicate_df, tax_df, treatment_target='catechin', day_min=7, day_max=35, top_n=3):
    """
    Melts, merges, and filters the expression data to return a dataframe 
    ready for plotting the top N active MAGs.
    """
    # 1. Melt to long format
    sample_cols = [col for col in df.columns if col.startswith('STM')]
    df_long = df.melt(
        id_vars=['MAG'], 
        value_vars=sample_cols,
        var_name='Sample', 
        value_name='geTMM'
    )
    
    # 2. Merge taxonomy and create the unique identifier string
    tax_subset = tax_df[['MAG', 'highest_host_tax_rank']].drop_duplicates()
    df_merged = df_long.merge(tax_subset, on='MAG', how='left')
    df_merged['unique_label'] = df_merged['MAG'].astype(str) + " | " + df_merged['highest_host_tax_rank'].astype(str)
    
    # 3. Merge metadata and filter for target treatment and day range
    df_merged = df_merged.merge(replicate_df[['Sample', 'treatment', 'day']], on='Sample', how='left')
    subset = df_merged[
        (df_merged['treatment'] == treatment_target) & 
        (df_merged['day'] >= day_min) & 
        (df_merged['day'] <= day_max)
    ].copy()
    
    # 4. Identify the top N MAGs based on total expression in this subset
    top_labels = subset.groupby('unique_label')['geTMM'].sum().nlargest(top_n).index
    
    # 5. Filter the dataframe to only include these top MAGs
    plot_data = subset[subset['unique_label'].isin(top_labels)].copy()
    
    return plot_data, top_labels

# --- 1. Prepare Data for Plotting ---

# Process Degraders (Assuming hydrogenase_uptake_catechin contains the degrader MAGs)
degraders_data, top_degrader_labels = get_top_mags_plot_data(
    hydrogenase_uptake_catechin, 
    replicate_frame, 
    MAG_metaT, 
    top_n=3
)

# Process Non-Degraders
non_degraders_data, top_non_degrader_labels = get_top_mags_plot_data(
    hydrogenase_uptake_non_degraders, 
    replicate_frame, 
    MAG_metaT, 
    top_n=3
)

# --- 2. Generate Box Plots ---

# Set up the matplotlib figure with two subplots (2 rows, 1 column)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))

# Custom styling for seaborn
sns.set_theme(style="whitegrid")

# Plot A: Top 3 Degraders
sns.boxplot(
    data=degraders_data,
    x='geTMM',
    y='unique_label',
    ax=axes[0],
    order=top_degrader_labels, # Ensures they are plotted from highest to lowest expression
    palette='viridis',
    showfliers=True
)
axes[0].set_title('Top 3 Degraders: Hydrogenase Uptake Activity (Catechin, Days 7-35)', fontsize=14, fontweight='bold')
axes[0].set_xlabel('Expression (geTMM)', fontsize=12)
axes[0].set_ylabel('MAG | Taxonomy', fontsize=12)

# Plot B: Top 3 Non-Degraders
sns.boxplot(
    data=non_degraders_data,
    x='geTMM',
    y='unique_label',
    ax=axes[1],
    order=top_non_degrader_labels,
    palette='magma',
    showfliers=True
)
axes[1].set_title('Top 3 Non-Degraders: Hydrogenase Uptake Activity (Catechin, Days 7-35)', fontsize=14, fontweight='bold')
axes[1].set_xlabel('Expression (geTMM)', fontsize=12)
axes[1].set_ylabel('MAG | Taxonomy', fontsize=12)

# Adjust layout to prevent overlapping labels
plt.tight_layout()

# Display the plot
plt.show()


# In[93]:


camper_genes_df.loc[camper_genes_df['highest_host_tax_rank'] == 'g__JAEXAI01'].camper_definition.unique()


# In[94]:


camper_genes_df.loc[camper_genes_df['highest_host_tax_rank'] == 'g__JAATFL01'].camper_definition.unique()


# In[95]:


camper_genes_df.columns


# # Plot abundance of all CAMPER pathways to see if Lignan degradation also increased

# In[96]:


import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy import stats

def get_significance_stars(p_val):
    """Returns stars based on p-value significance levels."""
    if p_val < 0.001:
        return "***"
    elif p_val < 0.01:
        return "**"
    elif p_val < 0.05:
        return "*"
    else:
        return "ns"

def plot_taxon_pathways_title_stats(df, taxon_name, replicate_frame):
    # 1. Filter for the specific taxon
    taxon_df = df[df['highest_host_tax_rank'] == taxon_name].copy()
    
    if taxon_df.empty:
        print(f"No data found for {taxon_name}")
        return

    # 2. Reshape from Wide to Long format
    sample_cols = [col for col in taxon_df.columns if col.startswith('STM')]
    df_long = taxon_df.melt(
        id_vars=['gene', 'camper_definition'], 
        value_vars=sample_cols,
        var_name='Sample', 
        value_name='geTMM'
    )

    # 3. Merge with replicate_frame
    df_plot = df_long.merge(replicate_frame[['Sample', 'treatment', 'day']], on='Sample', how='left')

    # 4. SUM PER SAMPLE
    df_summed = df_plot.groupby(['Sample', 'treatment', 'day', 'camper_definition'])['geTMM'].sum().reset_index()

    # 5. Identify unique definitions and grid size
    unique_defs = df_summed['camper_definition'].unique()
    num_defs = len(unique_defs)
    cols = 3
    rows = math.ceil(num_defs / cols)

    # 6. Plotting
    fig, axes = plt.subplots(rows, cols, figsize=(16, rows * 4), sharex=True)
    axes = axes.flatten()
    
    colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

    for i, gene_def in enumerate(unique_defs):
        ax = axes[i]
        gene_data = df_summed[df_summed['camper_definition'] == gene_def]
        
        # --- Statistical Test: Welch's t-test (Day 14-35) ---
        stats_window = gene_data[gene_data['day'] >= 14]
        cat_vals = stats_window[stats_window['treatment'] == 'catechin']['geTMM']
        una_vals = stats_window[stats_window['treatment'] == 'unamended']['geTMM']
        
        if len(cat_vals) > 1 and len(una_vals) > 1:
            _, p_val = stats.ttest_ind(cat_vals, una_vals, equal_var=False)
            sig_label = get_significance_stars(p_val)
        else:
            sig_label = "N/A"

        # Calculate mean for trend lines
        summary = gene_data.groupby(['treatment', 'day'])['geTMM'].agg(mean='mean').reset_index()

        for treatment in ['catechin', 'unamended']:
            t_data = gene_data[gene_data['treatment'] == treatment]
            t_summary = summary[summary['treatment'] == treatment]
            color = colors.get(treatment, 'gray')

            ax.scatter(t_data['day'], t_data['geTMM'], color=color, alpha=0.4, s=30, zorder=3)
            ax.plot(t_summary['day'], t_summary['mean'], color=color, linewidth=2.5, label=treatment, zorder=4)

        # --- Formatting with Significance in Title ---
        gene_short_name = gene_def.split(';')[0]
        ax.set_title(f"{gene_short_name} ({sig_label})", fontsize=12, fontweight='bold')
        
        ax.set_ylabel('Total geTMM', fontsize=9, fontweight='bold')
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.set_xticks([0, 7, 14, 21, 35])
        
        if i == 0: 
            ax.legend(title='Treatment', fontsize=9)
        if i >= (rows - 1) * cols: 
            ax.set_xlabel('Day', fontweight='bold')

    # Cleanup unused axes
    for j in range(i + 1, len(axes)): 
        axes[j].axis('off')

    plt.suptitle(f'Polyphenol Pathway Expression: {taxon_name}\nWelch\'s t-test (Days 14-35): *p<0.05, **p<0.01, ***p<0.001', 
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/{taxon_name}_polyphenol_pathway_expression.pdf', dpi=300)
    plt.show()

# --- Execute ---
plot_taxon_pathways_title_stats(camper_genes_df, 'g__JAEXAI01', replicate_frame)
plot_taxon_pathways_title_stats(camper_genes_df, 'g__JAATFL01', replicate_frame)


# In[97]:


import math
import numpy as np
import pandas as pd
from scipy import stats

def generate_s12_supplemental_table(df, taxa_list, replicate_frame, output_dir):
    """
    Processes multiple taxa data, extracts long-form expressions, 
    calculates window-restricted statistical tests, and saves a 
    clean, tidy supplemental CSV.
    """
    all_compiled_records = []
    
    for taxon_name in taxa_list:
        # 1. Filter for specific taxon
        taxon_df = df[df['highest_host_tax_rank'] == taxon_name].copy()
        if taxon_df.empty:
            print(f"Skipping {taxon_name}: No data found.")
            continue

        # 2. Reshape from Wide to Long format
        sample_cols = [col for col in taxon_df.columns if col.startswith('STM')]
        df_long = taxon_df.melt(
            id_vars=['gene', 'camper_definition'], 
            value_vars=sample_cols,
            var_name='Sample', 
            value_name='geTMM'
        )

        # 3. Merge with replicate metadata
        df_plot = df_long.merge(replicate_frame[['Sample', 'treatment', 'day']], on='Sample', how='left')

        # 4. Sum per individual sample profile
        df_summed = df_plot.groupby(['Sample', 'treatment', 'day', 'camper_definition'])['geTMM'].sum().reset_index()
        unique_defs = df_summed['camper_definition'].unique()

        # 5. Process each pathway definition and capture metrics
        for gene_def in unique_defs:
            gene_data = df_summed[df_summed['camper_definition'] == gene_def]
            
            # Statistical Test window constraints: Days 14-35 only
            stats_window = gene_data[gene_data['day'] >= 14]
            cat_vals = stats_window[stats_window['treatment'] == 'catechin']['geTMM']
            una_vals = stats_window[stats_window['treatment'] == 'unamended']['geTMM']
            
            p_val_numeric = np.nan
            if len(cat_vals) > 1 and len(una_vals) > 1:
                _, p_val = stats.ttest_ind(cat_vals, una_vals, equal_var=False)
                p_val_numeric = round(p_val, 6)

            # Build row dictionaries for each replicate point
            for _, row in gene_data.iterrows():
                gene_short_name = gene_def.split(';')[0]
                all_compiled_records.append({
                    'Taxon_Host_Rank': taxon_name,
                    'Camper_Definition_Full': gene_def,
                    'Gene_Short_Name': gene_short_name,
                    'Sample_ID': row['Sample'],
                    'Treatment': row['treatment'].capitalize(),
                    'Timepoint_Days': int(row['day']),
                    'Total_Expression_geTMM': round(row['geTMM'], 4),
                    'Welch_ttest_p_value_Days14_35': p_val_numeric
                })

    # 6. Convert to DataFrame and organize layout
    df_supplemental_s12 = pd.DataFrame(all_compiled_records)
    
    # Sort logically so reviewers can query by organism, then functional pathway matrix
    sorting_order = ['Taxon_Host_Rank', 'Gene_Short_Name', 'Treatment', 'Timepoint_Days', 'Sample_ID']
    df_supplemental_s12 = df_supplemental_s12.sort_values(by=sorting_order).reset_index(drop=True)
    
    # 7. Export file
    import os
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'S12_pathway_expression_stats.csv')
    df_supplemental_s12.to_csv(csv_path, index=False)
    
    print(f"Supplemental data table successfully saved to: {csv_path}")
    return df_supplemental_s12

# --- Execute compilation and save to target directory ---
target_taxa = ['g__JAEXAI01', 'g__JAATFL01']
destination_directory = '/users/PAS1573/riddell26/data'

df_s12_final = generate_s12_supplemental_table(
    df=camper_genes_df, 
    taxa_list=target_taxa, 
    replicate_frame=replicate_frame, 
    output_dir=destination_directory
)

df_s12_final


# # Check abundance of all these genes

# In[98]:


import matplotlib.pyplot as plt
import pandas as pd

# Define the directions (Selecting first 4 for a 2x2 grid)
hydro_directions = ['uptake', 'bifurcating', 'sensing', 'producing']
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Create the 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()  # Flatten to iterate easily

for i, direction in enumerate(hydro_directions):
    ax = axes[i]
    
    # 1. Filter using df_merged
    dir_only_df = df_merged.loc[ 
        (df_merged['enzyme'] == direction)
    ].copy()
    
    if dir_only_df.empty:
        ax.set_title(f'{direction.capitalize()} (No Data)')
        continue

    # 2. Aggregation and Merge
    dir_total_per_sample = dir_only_df.groupby(['sample']).agg({'geTMM': 'sum'}).reset_index()
    dir_total_per_sample.columns = ['Sample', 'total_getmm']
    dir_total_per_sample = dir_total_per_sample.merge(replicate_frame, on='Sample', how='left')

    # Determine time column name
    time_col = 'day' if 'day' in dir_total_per_sample.columns else 'time'
    
    # 3. Calculate mean stats
    summary = dir_total_per_sample.groupby(['treatment', time_col])['total_getmm'].agg([
        ('mean', 'mean')
    ]).reset_index()

    # 4. Plotting
    for treatment in summary['treatment'].unique():
        t_summary = summary[summary['treatment'] == treatment]
        replicates = dir_total_per_sample[dir_total_per_sample['treatment'] == treatment]
        
        color = colors.get(treatment, 'gray')
        
        # Plot replicates
        ax.scatter(replicates[time_col], replicates['total_getmm'],
                   color=color, alpha=0.5, s=30, zorder=3)
        
        # Plot mean line
        ax.plot(t_summary[time_col], t_summary['mean'], 
                linewidth=2, label=treatment, color=color, zorder=4)
        
    # 5. Facet Formatting
    ax.set_title(f'Hydrogenase: {direction.capitalize()}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_xticks([0, 7, 14, 21, 35])
    
    # Only add labels to the outer edges
    if i >= 2: ax.set_xlabel('Day', fontweight='bold')
    if i % 2 == 0: ax.set_ylabel('Total geTMM', fontweight='bold')
    
    # Add legend only to the first plot to keep it clean
    if i == 0:
        ax.legend(title='Treatment', fontsize=9)

plt.tight_layout()

# Save the facet grid
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10B-Hydrogenase_Facet_Grid.pdf', dpi=300)
plt.show()


# In[99]:


sorted_mags


# In[100]:


all_enriched_set = enriched_3plus_timepoints.union(enriched_2_timepoints).union(enriched_1_timepoint)


# In[101]:


all_enriched_set - set(sorted_mags)


# In[102]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Define enzymes to include
enzymes = ['CLD', # Lignans
            'EDL', # Lignans
            'PhcD', # Lignans
            'BER', # Lignans
            'SesA', # Lignans
            'PinZ', # Lignans
            'PAH', # Lignans
            'GLM', # Lignans
            'PhcC', # Lignans
            'PhcG', # Lignans
            
            'EhyB', # Phenyl-propenes
            'CalA', # Phenyl-propenes
            'CalB', # Phenyl-propenes
            'IEM', # Phenyl-propenes
            
            'FdeE', # Flavonones
            'FdeC', # Flavonones
           
            'CurA', # Curcuminoids
            'GRAD', # Phenolic acids
            'sam5', # Isoflavonoids, cinnamic acids

            'ChlE', # Cinnamic acid metabolism

            'FCS', # Non-specific
            ]

# Define timepoints to include
timepoints = [14, 21, 35]

# Pseudocount to avoid divide-by-zero
pseudocount = 0.0001

def calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus, enriched_2, enriched_1, mag_to_tax, df_merged):
    """
    Calculate log2 fold-change (catechin/unamended) for each MAG, enzyme, and timepoint.
    Only includes MAGs enriched in at least one timepoint.
    
    Returns:
    --------
    fc_df : DataFrame
        Dataframe with columns for each enzyme-timepoint combination
    """
    
    # Combine all enriched MAGs and sort
    sorted_mags = sorted(enriched_3plus) + sorted(enriched_2) + sorted(enriched_1)
    
    # Create column names for each enzyme-timepoint combination
    column_names = []
    for enzyme in enzymes:
        for time in timepoints:
            column_names.append(f"{enzyme}_day{time}")
    
    # Initialize matrix to store fold-changes
    fc_matrix = np.full((len(sorted_mags), len(column_names)), np.nan)
    
    # Calculate fold-change for each MAG, enzyme, and timepoint
    for i, mag in enumerate(sorted_mags):
        col_idx = 0
        for enzyme in enzymes:
            for time in timepoints:
                # Get data for this MAG, enzyme, and timepoint
                mag_enzyme_time_data = df_merged[(df_merged['MAG'] == mag) & 
                                              (df_merged['enzyme'] == enzyme) & 
                                              (df_merged['time'] == time)]
                
                if len(mag_enzyme_time_data) == 0:
                    # No data - leave as NaN (will be black in heatmap)
                    col_idx += 1
                    continue
                
                # Calculate mean expression for each treatment at this timepoint
                catechin_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'catechin']['geTMM'].mean()
                unamended_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'unamended']['geTMM'].mean()
                
                # Handle NaN values
                if pd.isna(catechin_mean):
                    catechin_mean = 0
                if pd.isna(unamended_mean):
                    unamended_mean = 0
                
                # Calculate log2 fold-change with pseudocount
                log2_fc = np.log2((catechin_mean + pseudocount) / (unamended_mean + pseudocount))
                fc_matrix[i, col_idx] = log2_fc
                
                col_idx += 1
    
    # Create DataFrame
    fc_df = pd.DataFrame(fc_matrix, index=sorted_mags, columns=column_names)
    
    # Add taxonomic information
    fc_df['Taxonomy'] = sorted_mags
    
    # Add enrichment type
    enrichment_type = []
    for mag in sorted_mags:
        if mag in enriched_3plus:
            enrichment_type.append('3 timepoints')
        elif mag in enriched_2:
            enrichment_type.append('2 timepoints')
        elif mag in enriched_1:
            enrichment_type.append('1 timepoint')
    fc_df['Enrichment_Type'] = enrichment_type
    
    return fc_df

def plot_foldchange_heatmap_enriched_only(fc_df, mag_to_tax):
    """
    Create heatmap of log2 fold-changes with timepoints for each enzyme.
    Enzymes on y-axis, MAGs on x-axis, with timepoints as sub-columns.
    Only shows enriched MAGs.
    """
    
    # Separate annotation columns from data
    data_cols = [col for col in fc_df.columns if col not in ['Taxonomy', 'Enrichment_Type']]
    heatmap_data = fc_df[data_cols]
    
    # Transpose the data (enzymes will be rows)
    heatmap_data = heatmap_data.T
    
    # Split data by enrichment type for column organization
    enriched_3plus_idx = fc_df[fc_df['Enrichment_Type'] == '3 timepoints'].index
    enriched_2_idx = fc_df[fc_df['Enrichment_Type'] == '2 timepoints'].index
    enriched_1_idx = fc_df[fc_df['Enrichment_Type'] == '1 timepoint'].index
    
    # Reorder columns by enrichment type
    column_order = list(enriched_3plus_idx) + list(enriched_2_idx) + list(enriched_1_idx)
    heatmap_data = heatmap_data[column_order]
    
    # Get taxonomy labels and enrichment types in same order
    taxonomy_labels = [f"{mag}\n{mag_to_tax.get(mag, 'Unknown')}; {completeness_contamination_dict.get(mag, 'Unknown')}" for mag in column_order]
    enrichment_types = fc_df.loc[column_order, 'Enrichment_Type'].values
    
    # Create figure
    fig = plt.figure(figsize=(max(16, len(column_order) * 0.5), 10))
    ax = fig.add_subplot(111)
    
    # Define color map - use a colormap that doesn't include white
    # RdBu_r goes from red (positive) through light colors to blue (negative)
    # We'll use a different palette that avoids white in the middle
    from matplotlib.colors import LinearSegmentedColormap
    
    # Create custom colormap: dark blue -> light blue -> light red -> dark red
    # This avoids pure white in the data range
    colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
              '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
    cmap.set_bad(color='#e0e0e0')  # Light gray for missing data
    
    # Determine vmin and vmax for symmetric color scale
    vmax = np.nanmax(np.abs(heatmap_data.values))
    vmin = -vmax
    
    # Create custom row labels (just enzyme names, timepoints will be in the cells)
    row_labels = []
    for col in heatmap_data.index:
        enzyme = col.split('_')[0]
        day = col.split('_day')[1]
        row_labels.append(f"{enzyme}_{day}")
    
    # Plot heatmap
    sns.heatmap(heatmap_data, 
                cmap=cmap,
                center=0,
                vmin=vmin,
                vmax=vmax,
                ax=ax,
                cbar=False,
                yticklabels=row_labels,
                xticklabels=taxonomy_labels,
                linewidths=0.5,
                linecolor='lightgray')
    
    # Add thicker black lines to separate enzymes (every 4 rows = 4 timepoints)
    for i in range(len(enzymes)):
        if i > 0:  # Don't draw line before first enzyme
            ax.axhline(y=i*len(timepoints), color='black', linewidth=3)
    
    # Add thicker vertical lines to separate MAGs
    for i in range(len(column_order)):
        if i > 0:
            ax.axvline(x=i, color='black', linewidth=1.5)
    
    ax.set_ylabel('Enzyme_Day', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9, length=0)
    
    # Set x-axis labels (taxonomy) at top and rotate 45 degrees
    ax.tick_params(axis='x', labelsize=7, rotation=45, labelbottom=False, labeltop=True, 
                  bottom=False, top=True)
    ax.xaxis.set_label_position('top')
    
    # Adjust label alignment for readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    
    # Add white vertical lines to separate enrichment groups
    boundaries = []
    prev_type = None
    for i, enrich_type in enumerate(enrichment_types):
        if enrich_type != prev_type and prev_type is not None:
            boundaries.append(i)
        prev_type = enrich_type
    
    for boundary in boundaries:
        ax.axvline(x=boundary, color='white', linewidth=4, linestyle='-')
    
    # Add text labels for each group at the top
    group_starts = [0] + boundaries
    group_ends = boundaries + [len(column_order)]
    
    for start, end in zip(group_starts, group_ends):
        midpoint = (start + end) / 2
        enrich_type = enrichment_types[start]
        ax.text(midpoint, len(heatmap_data) + 1.2, f'Catechin-enriched\n{enrich_type}',
               va='bottom', ha='center', fontsize=10, fontweight='bold',
               rotation=0)
    
    # Add colorbar at bottom
    cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.015])
    
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('log2-fold-change catechin / unamended mean GeTMM', fontsize=11, fontweight='bold')
    cbar.ax.tick_params(labelsize=9)
    
    # plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_enzyme_response_heatmap_enriched_only.pdf',
    #             dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig

# Generate the heatmap with only enriched MAGs
fc_df_enriched = calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus_timepoints, 
                                                                  enriched_2_timepoints, 
                                                                  enriched_1_timepoint,
                                                                  mag_to_tax, df_merged)
fig = plot_foldchange_heatmap_enriched_only(fc_df_enriched, mag_to_tax)

# Display summary statistics
print(f"\nSummary Statistics (Enriched MAGs Only):")
print(f"Total MAGs shown: {len(fc_df_enriched)}")
print(f"Enriched 3+ timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '3 timepoints')}")
print(f"Enriched 2 timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '2 timepoints')}")
print(f"Enriched 1 timepoint: {sum(fc_df_enriched['Enrichment_Type'] == '1 timepoint')}")
print(f"\nLog2 FC range: {np.nanmin(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f} to {np.nanmax(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f}")


# In[106]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Define enzymes to include
enzymes = ['FLR', 'CHI', 'FCR', 'PHY', 'PGR', 'pgthAB', # Flavonoids
           
           'fefe', 'nife', # Hydrogenases
    
    
            'CLD', # Lignans
            'EDL', # Lignans
            'PhcD', # Lignans
            'BER', # Lignans
            'SesA', # Lignans
            'PinZ', # Lignans
            'PAH', # Lignans
            'GLM', # Lignans
            'PhcC', # Lignans
            'PhcG', # Lignans
            
            'EhyB', # Phenyl-propenes
            'CalA', # Phenyl-propenes
            'CalB', # Phenyl-propenes
            'IEM', # Phenyl-propenes
            
            'FdeE', # Flavonones
            'FdeC', # Flavonones
           
            'CurA', # Curcuminoids
            'carABCDE', 'mhpABCD', 'hpaAB', 'hpaDEFGH', 'GRAD', # Phenolic acids
            'sam5', # Isoflavonoids, cinnamic acids

            'ChlE', # Cinnamic acid metabolism

            'FCS', # Non-specific
            ]

# Define timepoints to include
timepoints = [14, 21, 35]

# Pseudocount to avoid divide-by-zero
pseudocount = 0.0001

def calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus, enriched_2, enriched_1, mag_to_tax, df_merged):
    """
    Calculate log2 fold-change (catechin/unamended) for each MAG, enzyme, and timepoint.
    Only includes MAGs enriched in at least one timepoint.
    
    Returns:
    --------
    fc_df : DataFrame
        Dataframe with columns for each enzyme-timepoint combination
    """
    
    # Combine all enriched MAGs and sort
    sorted_mags = sorted(enriched_3plus) + sorted(enriched_2) + sorted(enriched_1)
    
    # Create column names for each enzyme-timepoint combination
    column_names = []
    for enzyme in enzymes:
        for time in timepoints:
            column_names.append(f"{enzyme}_day{time}")
    
    # Initialize matrix to store fold-changes
    fc_matrix = np.full((len(sorted_mags), len(column_names)), np.nan)
    
    # Calculate fold-change for each MAG, enzyme, and timepoint
    for i, mag in enumerate(sorted_mags):
        col_idx = 0
        for enzyme in enzymes:
            for time in timepoints:
                # Get data for this MAG, enzyme, and timepoint
                mag_enzyme_time_data = df_merged[(df_merged['MAG'] == mag) & 
                                              (df_merged['enzyme'] == enzyme) & 
                                              (df_merged['time'] == time)]
                
                if len(mag_enzyme_time_data) == 0:
                    # No data - leave as NaN (will be black in heatmap)
                    col_idx += 1
                    continue
                
                # Calculate mean expression for each treatment at this timepoint
                catechin_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'catechin']['geTMM'].mean()
                unamended_mean = mag_enzyme_time_data[mag_enzyme_time_data['treatment'] == 'unamended']['geTMM'].mean()
                
                # Handle NaN values
                if pd.isna(catechin_mean):
                    catechin_mean = 0
                if pd.isna(unamended_mean):
                    unamended_mean = 0
                
                # Calculate log2 fold-change with pseudocount
                log2_fc = np.log2((catechin_mean + pseudocount) / (unamended_mean + pseudocount))
                fc_matrix[i, col_idx] = log2_fc
                
                col_idx += 1
    
    # Create DataFrame
    fc_df = pd.DataFrame(fc_matrix, index=sorted_mags, columns=column_names)
    
    # Add taxonomic information
    fc_df['Taxonomy'] = sorted_mags
    
    # Add enrichment type
    enrichment_type = []
    for mag in sorted_mags:
        if mag in enriched_3plus:
            enrichment_type.append('3 timepoints')
        elif mag in enriched_2:
            enrichment_type.append('2 timepoints')
        elif mag in enriched_1:
            enrichment_type.append('1 timepoint')
    fc_df['Enrichment_Type'] = enrichment_type
    
    return fc_df

def plot_foldchange_heatmap_enriched_only(fc_df, mag_to_tax):
    """
    Create heatmap of log2 fold-changes with timepoints for each enzyme.
    Enzymes on y-axis, MAGs on x-axis, with timepoints as sub-columns.
    Only shows enriched MAGs.
    """
    
    # Separate annotation columns from data
    data_cols = [col for col in fc_df.columns if col not in ['Taxonomy', 'Enrichment_Type']]
    heatmap_data = fc_df[data_cols]
    
    # Transpose the data (enzymes will be rows)
    heatmap_data = heatmap_data.T
    
    # Split data by enrichment type for column organization
    enriched_3plus_idx = fc_df[fc_df['Enrichment_Type'] == '3 timepoints'].index
    enriched_2_idx = fc_df[fc_df['Enrichment_Type'] == '2 timepoints'].index
    enriched_1_idx = fc_df[fc_df['Enrichment_Type'] == '1 timepoint'].index
    
    # Reorder columns by enrichment type
    column_order = list(enriched_3plus_idx) + list(enriched_2_idx) + list(enriched_1_idx)
    heatmap_data = heatmap_data[column_order]
    
    # Get taxonomy labels and enrichment types in same order
    taxonomy_labels = [f"{mag}\n{mag_to_tax.get(mag, 'Unknown')}; {completeness_contamination_dict.get(mag, 'Unknown')}" for mag in column_order]
    enrichment_types = fc_df.loc[column_order, 'Enrichment_Type'].values
    
    # Create figure
    fig = plt.figure(figsize=(max(16, len(column_order) * 0.5), 20))
    ax = fig.add_subplot(111)
    
    # Define color map - use a colormap that doesn't include white
    # RdBu_r goes from red (positive) through light colors to blue (negative)
    # We'll use a different palette that avoids white in the middle
    from matplotlib.colors import LinearSegmentedColormap
    
    # Create custom colormap: dark blue -> light blue -> light red -> dark red
    # This avoids pure white in the data range
    colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
              '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
    cmap.set_bad(color='#e0e0e0')  # Light gray for missing data
    
    # Determine vmin and vmax for symmetric color scale
    vmax = np.nanmax(np.abs(heatmap_data.values))
    vmin = -vmax
    
    # Create custom row labels (just enzyme names, timepoints will be in the cells)
    row_labels = []
    for col in heatmap_data.index:
        enzyme = col.split('_')[0]
        day = col.split('_day')[1]
        row_labels.append(f"{enzyme}_{day}")
    
    # Plot heatmap
    sns.heatmap(heatmap_data, 
                cmap=cmap,
                center=0,
                vmin=vmin,
                vmax=vmax,
                ax=ax,
                cbar=False,
                yticklabels=row_labels,
                xticklabels=taxonomy_labels,
                linewidths=0.5,
                linecolor='lightgray')
    
    # Add thicker black lines to separate enzymes (every 4 rows = 4 timepoints)
    for i in range(len(enzymes)):
        if i > 0:  # Don't draw line before first enzyme
            ax.axhline(y=i*len(timepoints), color='black', linewidth=3)
    
    # Add thicker vertical lines to separate MAGs
    for i in range(len(column_order)):
        if i > 0:
            ax.axvline(x=i, color='black', linewidth=1.5)
    
    ax.set_ylabel('Enzyme_Day', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9, length=0)
    
    # Set x-axis labels (taxonomy) at top and rotate 45 degrees
    ax.tick_params(axis='x', labelsize=7, rotation=45, labelbottom=False, labeltop=True, 
                  bottom=False, top=True)
    ax.xaxis.set_label_position('top')
    
    # Adjust label alignment for readability
    plt.setp(ax.get_xticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    
    # Add white vertical lines to separate enrichment groups
    boundaries = []
    prev_type = None
    for i, enrich_type in enumerate(enrichment_types):
        if enrich_type != prev_type and prev_type is not None:
            boundaries.append(i)
        prev_type = enrich_type
    
    for boundary in boundaries:
        ax.axvline(x=boundary, color='white', linewidth=4, linestyle='-')
    
    # Add text labels for each group at the top
    group_starts = [0] + boundaries
    group_ends = boundaries + [len(column_order)]
    
    for start, end in zip(group_starts, group_ends):
        midpoint = (start + end) / 2
        enrich_type = enrichment_types[start]
        ax.text(midpoint, len(heatmap_data) + 1.2, f'Catechin-enriched\n{enrich_type}',
               va='bottom', ha='center', fontsize=10, fontweight='bold',
               rotation=0)
    
    # Add colorbar at bottom
    cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.015])
    
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('log2-fold-change catechin / unamended mean GeTMM', fontsize=11, fontweight='bold')
    cbar.ax.tick_params(labelsize=9)
    
    plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S10_enzyme_response_heatmap_enriched_only_additional_metabolisms.pdf',
                dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig

# Generate the heatmap with only enriched MAGs
fc_df_enriched = calculate_log2_foldchange_heatmap_enriched_only(enriched_3plus_timepoints, 
                                                                  enriched_2_timepoints, 
                                                                  enriched_1_timepoint,
                                                                  mag_to_tax, df_merged)
fig = plot_foldchange_heatmap_enriched_only(fc_df_enriched, mag_to_tax)

# Display summary statistics
print(f"\nSummary Statistics (Enriched MAGs Only):")
print(f"Total MAGs shown: {len(fc_df_enriched)}")
print(f"Enriched 3+ timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '3 timepoints')}")
print(f"Enriched 2 timepoints: {sum(fc_df_enriched['Enrichment_Type'] == '2 timepoints')}")
print(f"Enriched 1 timepoint: {sum(fc_df_enriched['Enrichment_Type'] == '1 timepoint')}")
print(f"\nLog2 FC range: {np.nanmin(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f} to {np.nanmax(fc_df_enriched[[col for col in fc_df_enriched.columns if col not in ['Taxonomy', 'Enrichment_Type']]].values):.2f}")


# In[107]:


fc_df_enriched.to_csv('/users/PAS1573/riddell26/data/S13_extended_metabolism_heatmap.csv')


# In[100]:


def analyze_additional_enzyme_enrichment(df, additional_enzymes):
    """
    Analyze enrichment of additional enzymes in catechin-amended samples for each MAG.
    
    Parameters:
    df: DataFrame with columns MAG, enzyme, treatment, time, geTMM, highest_host_tax_rank
    additional_enzymes: List of enzyme names to check
    
    Returns:
    Dictionary with MAG as key and dict of enriched enzymes and their enrichment counts
    """
    
    # Get all unique MAGs
    all_mags = df['MAG'].unique()
    
    # Create a mapping of MAG to taxonomic rank
    mag_to_tax = df.groupby('MAG')['highest_host_tax_rank'].first().to_dict()
    
    # Dictionary to store results for each MAG
    mag_results = {}
    
    # Analyze each MAG
    for mag in all_mags:
        enriched_enzymes_dict = {}
        
        # Check each additional enzyme
        for enzyme in additional_enzymes:
            enriched_timepoints = 0
            
            # Get data for this MAG and enzyme
            mag_enzyme_data = df[(df['MAG'] == mag) & (df['enzyme'] == enzyme)]
            
            # Skip if no data for this enzyme
            if len(mag_enzyme_data) == 0:
                continue
            
            # Check days 14, 21, 35
            for time in [14, 21, 35]:
                catechin_vals = mag_enzyme_data[(mag_enzyme_data['treatment'] == 'catechin') & 
                                               (mag_enzyme_data['time'] == time)]['geTMM'].values
                unamended_vals = mag_enzyme_data[(mag_enzyme_data['treatment'] == 'unamended') & 
                                                (mag_enzyme_data['time'] == time)]['geTMM'].values
                
                # Check if we have data for both treatments
                if len(catechin_vals) > 0 and len(unamended_vals) > 0:
                    max_catechin = catechin_vals.max()
                    mean_catechin = catechin_vals.mean()
                    max_unamended = unamended_vals.max()
                    mean_unamended = unamended_vals.mean()
                    
                    # Check for enrichment: BOTH max catechin > max unamended AND mean catechin > mean unamended
                    if max_catechin > max_unamended and mean_catechin > mean_unamended:
                        enriched_timepoints += 1
            
            # If this enzyme was enriched at any timepoint, record it
            if enriched_timepoints > 0:
                enriched_enzymes_dict[enzyme] = enriched_timepoints
        
        # Store results for this MAG
        mag_results[mag] = {
            'taxonomy': mag_to_tax.get(mag, 'Unknown'),
            'enriched_enzymes': enriched_enzymes_dict,
            'num_enriched_enzymes': len(enriched_enzymes_dict)
        }
    
    return mag_results

df_enriched = df_merged.loc[df_merged.MAG.isin(enriched_mags)]

# Run the analysis
results = analyze_additional_enzyme_enrichment(df_enriched, additional_enzymes)

# Print results for each MAG
print("="*80)
print("ENRICHMENT ANALYSIS: Additional Enzymes in Catechin-Amended Samples")
print("="*80)
print()

# Count MAGs with at least one enriched enzyme
mags_with_enrichment = 0

# Sort MAGs by number of enriched enzymes (descending), then by MAG name
sorted_mags = sorted(results.keys(), 
                     key=lambda x: (results[x]['num_enriched_enzymes'], x), 
                     reverse=True)

for mag in sorted_mags:
    result = results[mag]
    
    if result['num_enriched_enzymes'] > 0:
        mags_with_enrichment += 1
        
    print(f"MAG: {mag}")
    print(f"  Taxonomy: {result['taxonomy']}")
    print(f"  Number of enriched enzymes: {result['num_enriched_enzymes']}")
    
    if result['enriched_enzymes']:
        print(f"  Enriched enzymes:")
        # Sort enzymes by enrichment count (descending)
        sorted_enzymes = sorted(result['enriched_enzymes'].items(), 
                               key=lambda x: x[1], 
                               reverse=True)
        for enzyme, count in sorted_enzymes:
            print(f"    - {enzyme}: enriched at {count} timepoint(s)")
    else:
        print(f"  No enriched enzymes found")
    print()

# Summary statistics
print("="*80)
print("SUMMARY")
print("="*80)
print(f"Total MAGs analyzed: {len(results)}")
print(f"MAGs with at least one enriched enzyme: {mags_with_enrichment}")
print(f"MAGs with no enriched enzymes: {len(results) - mags_with_enrichment}")

# Count unique taxonomies with enrichment
unique_taxa_with_enrichment = set()
for mag, result in results.items():
    if result['num_enriched_enzymes'] > 0:
        unique_taxa_with_enrichment.add(result['taxonomy'])

print(f"\nUnique taxonomic ranks with at least one enriched enzyme: {len(unique_taxa_with_enrichment)}")
print(f"Taxa: {sorted(unique_taxa_with_enrichment)}")


# In[1]:


get_ipython().system('jupyter nbconvert --to script 05-AB_S11-13_polyphenol-expression.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/05-AB_S11-13_polyphenol-expression')


# In[ ]:





# In[ ]:




