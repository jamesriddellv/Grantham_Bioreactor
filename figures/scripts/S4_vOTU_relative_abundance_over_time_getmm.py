#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[4]:


active_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/active_vOTUs_1_read_mapped_relative_abundance.tsv', sep='\t')
active_df.head()


# In[5]:


# Add gene key so we can filter the htseq table
gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/prodigal-gv/combined_manual_filtered_gene_ids.txt', sep='\t')
gff_df = gff_df[['vOTU', 'gene']]
gff_df['gene_position'] = gff_df['gene'].apply(lambda x: x.split('_')[1])
gff_df.head()


# In[6]:


# Read readcount data
data = pd.read_csv(
    '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED_no0s.tsv',
    sep='\t', header=None
)

# Assign column names
columns = ["gene",
  "STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004",
  "STM_0716_E_M_E025", "STM_0716_E_M_E027", "STM_0716_E_M_E029",
  "STM_0716_E_M_E030", "STM_0716_E_M_E031", "STM_0716_E_M_E033",
  "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
  "STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E054",
  "STM_0716_E_M_E055", "STM_0716_E_M_E056", "STM_0716_E_M_E058",
  "STM_0716_E_M_E059", "STM_0716_E_M_E060", "STM_0716_E_M_E062",
  "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E066",
  "STM_0716_E_M_E067", "STM_0716_E_M_E068", "STM_0716_E_M_E070",
  "STM_0716_E_M_E071", "STM_0716_E_M_E072", "STM_0716_E_M_E121",
  "STM_0716_E_M_E122", "STM_0716_E_M_E123", "STM_0716_E_M_E125",
  "STM_0716_E_M_E126", "STM_0716_E_M_E127", "STM_0716_E_M_E129",
  "STM_0716_E_M_E130", "STM_0716_E_M_E131"]
data.columns = columns

# Drop condensed tannin
data = data.drop(columns=['STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                       'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                       'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                       'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127'])

# Drop last two metadata rows
data = data.iloc[:-2, :]

# Add vOTU
data = data.merge(gff_df, on='gene', how='left')

# filter for only active vOTUs
data = data.loc[data['vOTU'].isin(active_df['vOTU'].unique())]



# In[7]:


vOTU_readcount = data.drop(columns=['gene', 'gene_position']).groupby('vOTU').sum()
vOTU_readcount


# In[8]:


df = vOTU_readcount


# # Species accumulation curve

# In[54]:


from itertools import permutations

# === Load your data ===
# === Transpose so columns are species ===
df_T = df.T  # Now rows = samples, columns = species

# === Updated species_accumulation function to return all data ===
def species_accumulation(df, permutations_count=100):
    species_acc = np.zeros((permutations_count, len(df)))

    for i in range(permutations_count):
        order = np.random.permutation(df.index)
        seen_species = set()
        counts = []

        for j, sample in enumerate(order):
            sample_species = df.loc[sample]
            new_species = sample_species[sample_species > 0].index
            seen_species.update(new_species)
            counts.append(len(seen_species))
        
        species_acc[i, :] = counts

    mean_acc = species_acc.mean(axis=0)
    min_acc = species_acc.min(axis=0)
    max_acc = species_acc.max(axis=0)
    
    return mean_acc, min_acc, max_acc

# === Compute curve data ===
mean_curve, min_curve, max_curve = species_accumulation(df_T)

# === Plot with min-max shaded range ===
sns.set(style="whitegrid")
x = np.arange(1, len(mean_curve) + 1)

ax = plt.figure(figsize=(6, 6))
sns.lineplot(x=x, y=mean_curve, label="Mean species richness", color='black')
plt.fill_between(x, min_curve, max_curve, alpha=0.3, label="Min–Max range", color='grey')

# Axis labels with larger font
plt.xlabel("Number of samples", fontsize=14)
plt.ylabel("vOTUs detected (cumulative)", fontsize=14)
xtick_positions = np.arange(1, 27, 5)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(fontsize=14)

# Start x-axis at 1
plt.xlim(left=1)

# Move legend to the right
# plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0, fontsize=12)
plt.legend().remove()
plt.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/01-B_specaccum_active_1_read_mapped.pdf", dpi=300)
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/01-B_specaccum_active_1_read_mapped.png", dpi=300)
plt.show()




# # Make for Catechin and Unamended

# In[55]:


# Split by catechin and unamended samples
replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[56]:


unamended_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'unamended') | (replicate_frame.day == 0)]['Sample'])
catechin_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'catechin') | (replicate_frame.day == 0)]['Sample'])


# In[57]:


df_T_U = df_T.loc[df_T.index.isin(unamended_samples)]
df_T_C = df_T.loc[df_T.index.isin(catechin_samples)]


# In[58]:


from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === Updated species_accumulation function to return all data ===
def species_accumulation(df, permutations_count=100):
    species_acc = np.zeros((permutations_count, len(df)))
    for i in range(permutations_count):
        order = np.random.permutation(df.index)
        seen_species = set()
        counts = []
        for j, sample in enumerate(order):
            sample_species = df.loc[sample]
            new_species = sample_species[sample_species > 0].index
            seen_species.update(new_species)
            counts.append(len(seen_species))
        
        species_acc[i, :] = counts
    mean_acc = species_acc.mean(axis=0)
    min_acc = species_acc.min(axis=0)
    max_acc = species_acc.max(axis=0)
    
    return mean_acc, min_acc, max_acc

# === Compute curve data for both treatments ===
mean_curve_U, min_curve_U, max_curve_U = species_accumulation(df_T_U)
mean_curve_C, min_curve_C, max_curve_C = species_accumulation(df_T_C)

# === Plot with min-max shaded range for both treatments ===
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))

# Unamended treatment
x_U = np.arange(1, len(mean_curve_U) + 1)
sns.lineplot(x=x_U, y=mean_curve_U, label="Unamended", color='#7bccc4', linewidth=2)
plt.fill_between(x_U, min_curve_U, max_curve_U, alpha=0.3, color='#7bccc4')

# Catechin treatment
x_C = np.arange(1, len(mean_curve_C) + 1)
sns.lineplot(x=x_C, y=mean_curve_C, label="Catechin", color='#fe9929', linewidth=2)
plt.fill_between(x_C, min_curve_C, max_curve_C, alpha=0.3, color='#fe9929')

# Axis labels with larger font
plt.xlabel("Number of samples", fontsize=14)
plt.ylabel("vOTUs detected (cumulative)", fontsize=14)

# Set x-ticks based on the maximum number of samples
max_samples = max(len(df_T_U), len(df_T_C))
xtick_positions = np.arange(1, max_samples + 1, 5)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(fontsize=14)

# Start x-axis at 1
plt.xlim(left=1)

# Add legend
plt.legend(fontsize=12, loc='lower right')
plt.grid(alpha=0.2)
plt.tight_layout()

# Uncomment to save figures
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/S1X_specaccum_active_by_treatment_1_read_mapped.pdf", dpi=300)
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/S1X_specaccum_active_by_treatment_1_read_mapped.png", dpi=300)
plt.show()


# In[59]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.loc[iphop_df['MAG_is_active'] == True]


# In[60]:


# is all lachnoclostridium also predicted to infect clostridium?
clos_df = (
    iphop_df
    .loc[(iphop_df['vOTU'].isin(set(df.index.unique())))]
    .groupby(['vOTU'])
    .agg({'highest_host_tax_rank': ';'.join})
    .reset_index()

)
clos_df.loc[clos_df['highest_host_tax_rank'].str.contains('g__Lachnoclostridium')]


# In[61]:


iphop_df.loc[iphop_df['vOTU'].isin(['contig_1022629', 'contig_1093732'])]


# In[62]:


# due to changes in taxonomy, we will change g__Lachnoclostridium to g__Clostridium
# iphop_df['G'] = iphop_df['G'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
# iphop_df['highest_host_tax_rank'] = iphop_df['highest_host_tax_rank'].str.replace('g__Lachnoclostridium', 'g__Clostridium')

iphop_df = iphop_df.sort_values(by=['vOTU', 'Confidence score']).drop_duplicates('vOTU')
iphop_df = iphop_df.loc[iphop_df['vOTU'].isin(list(df.index))]
iphop_df.head()


# In[63]:


iphop_df.loc[iphop_df['vOTU'] == 'contig_591846']


# In[64]:


df = df.merge(iphop_df, on='vOTU', how='left')


# In[65]:


df.columns


# In[66]:


keep_cols = ['highest_host_tax_rank', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']


# In[67]:


df['highest_host_tax_rank'] = df['highest_host_tax_rank'].fillna('no_assignment')


# In[68]:


by_tax_rank = df.groupby('highest_host_tax_rank').sum().reset_index()[keep_cols]


# In[69]:


by_tax_rank_melted = by_tax_rank.melt(id_vars=['highest_host_tax_rank'], var_name='Sample', value_name='sum_mean_abundance')
by_tax_rank_melted


# In[70]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# In[71]:


by_tax_rank_melted = by_tax_rank_melted.merge(replicate_frame, on='Sample', how='left')


# In[72]:


# add reads per sample
reads_per_sample = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv', usecols=['Sample', 'total reads (R1+R2)'])
by_tax_rank_melted = by_tax_rank_melted.merge(reads_per_sample, on='Sample', how='left')


# In[73]:


by_tax_rank_melted


# In[74]:


# Add prop abundance
by_sample_sum_abundance = by_tax_rank_melted.groupby(['Sample']).agg({'sum_mean_abundance': 'sum'}).reset_index().rename(columns={'sum_mean_abundance': 'sample_sum'})


# In[75]:


by_tax_rank_melted = by_tax_rank_melted.merge(by_sample_sum_abundance, on='Sample', how='left')


# In[76]:


by_tax_rank_melted['prop_abundance'] = by_tax_rank_melted['sum_mean_abundance'] / by_tax_rank_melted['sample_sum']
by_tax_rank_melted['prop_total'] = by_tax_rank_melted['sum_mean_abundance'] / by_tax_rank_melted['total reads (R1+R2)']


# In[77]:


by_tax_rank_melted


# In[78]:


by_tax_rank_melted['day'] = by_tax_rank_melted['day'].astype(int)


# In[79]:


# Step 1: Make sure 'taxon' is the correct column name

# Step 2: Group by treatment_day_replicate and taxon, calculate mean prop_abundance
grouped = (
    by_tax_rank_melted.groupby(['treatment_day_replicate', 'highest_host_tax_rank'])['prop_abundance']
    .mean()
    .reset_index()
)

# Step 3: Find max prop_abundance per taxon across all groups
max_abundance = grouped.groupby('highest_host_tax_rank')['prop_abundance'].max()

# Step 4: Identify taxa to keep (those with at least one group ≥ 5%)
taxa_to_keep = max_abundance[max_abundance >= 0.05].index

print(taxa_to_keep)

# From Step 4 but added UBA-1794 because it was identified in >10% of reads.
# taxa_to_keep = ['g__Clostridium', 'g__Fen-1039', 'g__JAEXAI01', 'g__JAGFXR01',
#        'g__Paludibacter', 'g__Pseudomonas_E', 'g__Terracidiphilus',
#        'g__UBA1794', 'g__Cfx3-03', 'g__Enterobacter']

# taxa_to_keep merged with MAGs > 5%
# taxa_to_keep = ['no_assignment', 'g__Clostridium', 'g__JAGFXR01', 'g__Paludibacter', 'g__Pseudomonas_E', 'g__Terracidiphilus', 'g__Cfx3-03']

# Step 5: Create a new column in original df with taxon or "other"
by_tax_rank_melted['taxon_grouped'] = by_tax_rank_melted['highest_host_tax_rank'].where(by_tax_rank_melted['highest_host_tax_rank'].isin(taxa_to_keep), 'Other (<5%)')


# In[80]:


by_tax_rank_melted.columns


# In[81]:


by_tax_rank_melted_grouped = by_tax_rank_melted.groupby(['Sample', 'sum_mean_abundance',
       'sample_treatment_day', 'replicate', 'treatment', 'day',
       'treatment_day_replicate', 'prop_abundance', 'prop_total',
       'taxon_grouped']).agg(taxon_grouped_prop_abundance=('prop_abundance', 'sum'), taxon_grouped_prop_total=('prop_total', 'sum')).reset_index()


# In[82]:


# Step 1: Filter unamended day 0 rows
day0_unamended = by_tax_rank_melted_grouped[
    (by_tax_rank_melted_grouped['treatment'] == 'unamended') &
    (by_tax_rank_melted_grouped['day'] == 0)
].copy()

# Step 2: Modify treatment to 'catechin'
day0_unamended['treatment'] = 'catechin'

# Optional: also update identifiers like 'sample_treatment_day' or 'treatment_day_replicate' if needed
# For example:
# day0_unamended['sample_treatment_day'] = day0_unamended['sample_treatment_day'].str.replace('unamended', 'catechin')
# day0_unamended['treatment_day_replicate'] = day0_unamended['treatment_day_replicate'].str.replace('unamended', 'catechin')

# Step 3: Append to original dataframe
by_tax_rank_melted_grouped = pd.concat([by_tax_rank_melted_grouped, day0_unamended], ignore_index=True)


# In[83]:


by_tax_rank_melted_grouped


# In[84]:


mag_color_dict = {
    'Other (<10%)': 'lightgray',
    'g_Clostridium': '#FF6B35',      # Vibrant orange
    'g_Fen-1039': '#004E89',         # Deep blue
    'g_JAEXAI01': '#00C9A7',         # Turquoise
    'g_JAGFXR01': '#845EC2',         # Purple
    'g_Paludibacter': '#FF1744',     # Bright red
    'g_Pseudomonas_E': '#FFD700',    # Gold
    'g_Terracidiphilus': '#32CD32',  # Lime green
    'g_UBA1794': '#FF69B4',           # Hot pink
    'g__Cfx3-03': '#8B4513',        # Saddle brown
    'g__Enterobacter': '#00CED1',   # Dark turquoise
}


# In[85]:


mag_color_dict = {
    'Other (<5%)': 'gray',
    'g__Clostridium': '#FF6B35',      # Vibrant orange
    'g__JAEXAI01': '#00C9A7',         # Turquoise
    'g__JAGFXR01': '#845EC2',         # Purple
    'g__PALSA-129': '#004E89',        # Deep blue
    'g__Paludibacter': '#FF1744',     # Bright red
    'g__Pseudomonas_E': '#FFD700',    # Gold
    'g__Terracidiphilus': '#32CD32',  # Lime green
    'g__UBA1794': '#FF69B4',          # Hot pink
    'no_assignment': 'lightgray',       # Saddle brown
}


# In[86]:


mag_color_dict = {
    'Other': 'gray',
    'g__Cfx3-03': '#8B4513',            # Saddle brown
    'g__Clostridium': '#9E0B0F',        # Sangria red
    'g__Enterobacter': '#00CED1',       # Dark turquoise
    'g__Fen-1039': '#FF8C00',           # Dark orange
    'g__Fen-1137': '#4682B4',           # Steel blue
    'g__Fen-1308': '#ADFF2F',           # Green yellow
    'g__JAAYUD01': '#9932CC',           # Dark orchid
    'g__JAGFXR01': '#845EC2',           # Purple
    'g__JAEXAI01': '#00C9A7',           # Turquoise
    'g__Janthinobacterium': '#1E90FF',  # Dodger blue
    'g__LD21': '#FF6B35',        # Vibrant orange
    'g__PALSA-129': '#004E89',          # Deep blue
    'g__Paludibacter': '#FF1744',       # Bright red
    'g__Pseudomonas_E': '#FFD700',      # Gold
    'g__Terracidiphilus': '#32CD32',    # Lime green
    'g__UBA1794': '#FF69B4',            # Hot pink
    'no_assignment': 'lightgrey'
}


# In[87]:


mag_color_dict = {
    'Other': 'gray',
    'g__Bacteroides': '#4682B4',        # Steel blue (new assignment)
    'g__Clostridium': '#9E0B0F',        # Sangria red (preserved)
    'g__JAEXAI01': '#00C9A7',           # Turquoise (preserved)
    'g__JAGFXR01': '#845EC2',           # Purple (preserved)
    'g__Lachnoclostridium': '#FF8C00',  # Dark orange (reassigned)
    'g__Paludibacter': '#FF1744',       # Bright red (preserved)
    'g__Pseudomonas_E': '#FFD700',      # Gold (preserved)
    'g__Terracidiphilus': '#32CD32',    # Lime green (preserved)
    'g__UBA1794': '#FF69B4',            # Hot pink (preserved)
    'no_assignment': 'lightgrey'        # Light grey (preserved)
}


# In[88]:


# Pivot the dataframe: rows = day, columns = taxon, values = mean prop abundance
df_grouped = (
    by_tax_rank_melted_grouped
    .groupby(['treatment', 'day', 'taxon_grouped'])['taxon_grouped_prop_abundance']
    .sum()
    .reset_index()
).sort_values(by='taxon_grouped_prop_abundance', ascending=False)

# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = df_grouped[df_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
df_pivot = df_pivot.fillna(0).sort_index()
df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("proportion of active vOTU of total RNA", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
plt.ylim(0, 1)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/stacked_area_{treatment_to_plot}_vOTU_metaT_simplified.png",
    dpi=300
)
plt.show()


# In[89]:


# Pivot the dataframe: rows = day, columns = taxon, values = mean prop abundance
df_grouped = (
    by_tax_rank_melted_grouped
    .groupby(['treatment', 'day', 'taxon_grouped'])['taxon_grouped_prop_total']
    .sum()
    .reset_index()
).sort_values(by='taxon_grouped_prop_total', ascending=False)

# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = df_grouped[df_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_total')
df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("proportion of active vOTU of total RNA", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks([0.0, 0.006], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
# plt.ylim(0, 1)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/stacked_area_{treatment_to_plot}_vOTU_metaT_TOTAL_simplified.png",
    dpi=300
)
plt.show()


# In[90]:


# Pivot the dataframe: rows = day, columns = taxon, values = mean prop abundance
df_grouped = (
    by_tax_rank_melted_grouped
    .groupby(['treatment', 'day', 'taxon_grouped'])['taxon_grouped_prop_abundance']
    .sum()
    .reset_index()
).sort_values(by='taxon_grouped_prop_abundance', ascending=False)

# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = df_grouped[df_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
df_pivot = df_pivot.fillna(0).sort_index()
df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("proportion of active vOTU total RNA", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
plt.ylim(0, 1)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/stacked_area_{treatment_to_plot}_vOTU_metaT_simplified.png",
    dpi=300
)
plt.show()


# In[91]:


# Pivot the dataframe: rows = day, columns = taxon, values = mean prop abundance
df_grouped = (
    by_tax_rank_melted_grouped
    .groupby(['treatment', 'day', 'taxon_grouped'])['taxon_grouped_prop_total']
    .sum()
    .reset_index()
).sort_values(by='taxon_grouped_prop_total', ascending=False)

# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = df_grouped[df_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_total')
df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("proportion of active vOTU of total RNA", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks([0.0, 0.01, 0.02, 0.03], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
# plt.ylim(0, 1)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4/stacked_area_{treatment_to_plot}_vOTU_metaT_TOTAL_simplified.png",
    dpi=300
)
plt.show()


# In[45]:


df.loc[df['highest_host_tax_rank'] == 'g__JAGFXR01']


# In[83]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__JAGFXR01'] * 100


# In[38]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__Pseudomonas_E'] * 100


# In[39]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__Clostridium'] * 100


# In[38]:


df_pivot.sum(axis=1)


# # Plot individual phages from these >10%

# In[44]:


import matplotlib.colors as mcolors
treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

# Function to create a lighter color
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by mixing it with white.
    `amount` controls the mix ratio (0 = no change, 1 = white).
    """
    rgb = np.array(mcolors.to_rgb(color))  # Convert to RGB
    white = np.array([1, 1, 1])  # White color
    return mcolors.to_hex((1 - amount) * rgb + amount * white)  # Blend with white


def plot_vOTU_abundance_by_treatment(highest_host_tax_rank, df, with_individual=False, has_recombinase=False, is_active=False, with_replicate_points=True):
    
    if highest_host_tax_rank != []:
        df_subset = df.loc[
            (df['highest_host_tax_rank'].isin(highest_host_tax_rank)) & 
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()

    else:
        df_subset = df.loc[
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()   

    if has_recombinase:
        df_subset = df_subset.loc[df['is_provirus_aggregated'] == True]

    if is_active:
        df_subset = df_subset.loc[df['is_active'] == True]
    
    sns.set_style("white")
    treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

    plt.figure(figsize=(8, 6))

    if with_individual:
        for (votu, treatment), group in df_subset.groupby(['vOTU', 'treatment']):
            # light_color = lighten_color(treatment_dict[treatment], amount=0.5)
            sns.lineplot(
                data=group, x='day', y='log2_abundance',
                color=treatment_dict[treatment], linewidth=0.8, alpha=0.5, errorbar=None, legend=False
            )


    # Plot mean lines per treatment
    for treatment, color in treatment_dict.items():
        df_treat = df_subset[df_subset['treatment'] == treatment]
        sns.lineplot(
            data=df_treat, x='day', y='log2_abundance',
            color=color, linewidth=2.5, errorbar=('ci', 95), label=treatment
        )


    if with_replicate_points and 'replicate' in df_subset.columns:
        # Plot replicate means as points
        mean_by_replicate = (
            df_subset.groupby(['treatment', 'day', 'replicate'])['log2_abundance']
            .mean()
            .reset_index()
        )
        for treatment, color in treatment_dict.items():
            subset = mean_by_replicate[mean_by_replicate['treatment'] == treatment]
            sns.scatterplot(
                data=subset, x='day', y='log2_abundance',
                color=color, s=60, edgecolor='black', linewidth=0.5,
                alpha=0.7, label=None
            )
    
    plt.xlabel("Day", fontsize=12)
    plt.ylabel("log2(Abundance)", fontsize=12)
    # plt.title(f'vOTU Abundance Over Time ({highest_host_tax_rank})', fontsize=14)
    plt.legend(title='Treatment')
    plt.legend().remove()
    plt.tight_layout()
    plt.show()


# In[45]:


df


# In[46]:


df.columns


# In[47]:


df_melted = df[['vOTU', 'highest_host_tax_rank', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']].melt(id_vars=['vOTU', 'highest_host_tax_rank'], var_name='Sample', value_name='abundance')

df_melted = df_melted.merge(replicate_frame, on='Sample', how='left')
df_melted['log2_abundance'] = df_melted['abundance'].apply(lambda x: np.log2(x + 1))
df_melted.head()


# In[48]:


plot_vOTU_abundance_by_treatment(['g__JAGFXR01'], df_melted, with_individual=True)


# In[49]:


plot_vOTU_abundance_by_treatment(['g__Paludibacter'], df_melted, with_individual=True)


# # label by catechin-enriched or catechin-sensitive

# In[ ]:





# In[ ]:




