#!/usr/bin/env python
# coding: utf-8

# # Generates Figure 1B and alternate versions of the plot (scaled by total GeTMM instead of proportion of total)

# In[ ]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from print_versions import print_versions
print_versions(globals())


# In[ ]:


# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


# In[ ]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/geTMM_table.csv')
df.head()


# In[ ]:


# get from DRAM table, load only the gene column and fasta column, and count number of genes.
dram_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/reactorEMERGE_annotations.tsv', sep='\t', usecols=['fasta', 'gene_position'])
dram_df


# In[ ]:


num_genes = dram_df.fasta.value_counts().reset_index()
num_genes.columns=['MAG', 'num_genes']
num_genes


# In[ ]:


df_grouped = df.drop(columns=['gene']).groupby('fasta').sum().reset_index()
df_grouped.rename(columns={'fasta': 'MAG'}, inplace=True)
df_grouped = df_grouped.merge(num_genes, on='MAG', how='left')
df_grouped.head()


# In[ ]:


# Divide columns 1 to 26 (since slicing in iloc is exclusive of the end index)
df_grouped.iloc[:, 1:27] = df_grouped.iloc[:, 1:27].div(df_grouped['num_genes'], axis=0)
df_grouped = df_grouped.drop(columns=['num_genes'])
df_grouped


# In[ ]:


# join active MAGs
gtdb_df = pd.read_excel('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_abundance_table.xlsx', usecols=['MAG', 'GTDB'])

df_grouped = df_grouped.merge(gtdb_df, on='MAG', how='left')
df_grouped


# In[ ]:


df_grouped[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = df_grouped['GTDB'].str.split(';', expand=True)


# In[ ]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

df_grouped['highest_host_tax_rank'] = df_grouped.apply(get_highest_resolution, axis=1)


# In[ ]:


# MAG relative abundance with GTDB taxonomy features
df_grouped.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_GTDB_getmm.tsv', sep='\t', index=False)


# In[ ]:


by_tax_rank = df_grouped.groupby(['highest_host_tax_rank']).sum().reset_index()
by_tax_rank = by_tax_rank[['highest_host_tax_rank', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003',
       'STM_0716_E_M_E004', 'STM_0716_E_M_E025', 'STM_0716_E_M_E027',
       'STM_0716_E_M_E033', 'STM_0716_E_M_E034', 'STM_0716_E_M_E035',
       'STM_0716_E_M_E050', 'STM_0716_E_M_E051', 'STM_0716_E_M_E052',
       'STM_0716_E_M_E058', 'STM_0716_E_M_E059', 'STM_0716_E_M_E060',
       'STM_0716_E_M_E062', 'STM_0716_E_M_E063', 'STM_0716_E_M_E064',
       'STM_0716_E_M_E070', 'STM_0716_E_M_E071', 'STM_0716_E_M_E072',
       'STM_0716_E_M_E121', 'STM_0716_E_M_E122', 'STM_0716_E_M_E123',
       'STM_0716_E_M_E129', 'STM_0716_E_M_E130', 'STM_0716_E_M_E131']]
by_tax_rank


# In[ ]:


by_tax_rank_melted = by_tax_rank.melt(id_vars=['highest_host_tax_rank'], var_name='Sample', value_name='getmm')
by_tax_rank_melted


# In[ ]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# In[ ]:


by_tax_rank_melted = by_tax_rank_melted.merge(replicate_frame, on='Sample', how='left')


# In[ ]:


by_tax_rank_melted.sort_values(by='getmm', ascending=False)


# In[ ]:


# Take the mean across replicates, since there are fewer replicates of unamended day 14
by_tax_rank_melted = by_tax_rank_melted.groupby(['highest_host_tax_rank', 'treatment', 'day']).agg({'getmm': 'mean'}).reset_index()


# In[ ]:


# Add prop abundance
by_sample_sum_abundance = by_tax_rank_melted.groupby(['treatment', 'day']).agg({'getmm': 'sum'}).reset_index().rename(columns={'getmm': 'total_getmm'})


# In[ ]:


by_tax_rank_melted = by_tax_rank_melted.merge(by_sample_sum_abundance, on=['treatment', 'day'], how='left')


# In[ ]:


by_tax_rank_melted['prop_abundance'] = by_tax_rank_melted['getmm'] / by_tax_rank_melted['total_getmm']


# In[ ]:


by_tax_rank_melted


# In[ ]:


by_tax_rank_melted['day'] = by_tax_rank_melted['day'].astype(int)


# In[ ]:


# Step 3: Find max prop_abundance per taxon across all groups
max_abundance = by_tax_rank_melted.groupby('highest_host_tax_rank')['prop_abundance'].max()

# Step 4: Identify taxa to keep (those with at least one group ≥ 5%)
taxa_to_keep = max_abundance[max_abundance >= 0.05].index
taxa_to_keep


# In[ ]:


# Step 5: Create a new column in original df with taxon or "other"
by_tax_rank_melted['taxon_grouped'] = by_tax_rank_melted['highest_host_tax_rank'].where(by_tax_rank_melted['highest_host_tax_rank'].isin(taxa_to_keep), 'Other')


# In[ ]:


by_tax_rank_melted.columns


# In[ ]:


by_tax_rank_melted


# In[ ]:


by_tax_rank_melted_grouped = (
    by_tax_rank_melted
    .groupby(['treatment', 'day', 'taxon_grouped'])
    .agg(taxon_grouped_getmm=('getmm', 'sum'), taxon_grouped_prop_abundance=('prop_abundance', 'sum'))
    .reset_index()
)
by_tax_rank_melted_grouped


# In[ ]:


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


# In[ ]:


mag_color_dict = {
    'Other': 'black',
    
    # g__Clostridium group (red shades)
    'g__Clostridium': '#9E0B0F',                    # Sangria red
    'g__Clostridium__contig_507316': '#DC143C',     # Crimson
    'g__Clostridium__contig_588865': '#FF4500',     # Orange red
    'g__Clostridium__contig_755337': '#B22222',     # Fire brick red
    
    # g__Fen-1039 group (orange shades)
    'g__Fen-1039': '#FF8C00',                       # Dark orange

    # g__Fen-1308 group (blue shade)
    'g__Fen-1308': '#1E90FF',                       # Dodger blue
    
    # g__JAGFXR01 group (purple shades)
    'g__JAGFXR01': '#845EC2',                       # Purple
    'g__JAGFXR01__contig_591846': '#5C4188',        # Dark purple
    
    # g__Cfx3-03 group (green shades)
    'g__Cfx3-03': '#228B22',                        # Forest green
    
    # g__Paludibacter group (red-pink shades)
    'g__Paludibacter': '#FF1744',                   # Bright red
    'g__Paludibacter__contig_724521': '#E91E63',    # Pink-red
    
    # g__Pseudomonas_E group (gold shades)
    'g__Pseudomonas_E': '#FFD700',                  # Gold
    
    # g__Terracidiphilus group (green shades)
    'g__Terracidiphilus': '#32CD32',                # Lime green
    
    # g__UBA1794 group (pink shades)
    'g__UBA1794': '#FF69B4',                        # Hot pink
    
    # no_assignment group (gray shades)
    'no_assignment__contig_1022629': '#C0C0C0',     # Silver
    'no_assignment__contig_1093732': '#D3D3D3',     # Light gray
    'no_assignment__contig_709158': '#696969',      # Dim gray
    'no_assignment__contig_713418': '#778899',      # Light slate gray
    'no_assignment__contig_773239': '#708090',      # Slate gray
    'no_assignment__contig_776767': '#556B2F',      # Dark olive green
    'no_assignment__contig_861080': '#2F4F4F',      # Dark slate gray
    'no_assignment': '#808080',                     # Medium gray
}


# In[ ]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_getmm')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("MAG relative abundance by genus (GeTMM)", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
# plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
# plt.ylim(0, 1)

# handles, labels = plt.gca().get_legend_handles_labels()
# plt.legend(reversed(handles), reversed(labels), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_MAG_stacked_area_{treatment_to_plot}.pdf", dpi=300)
plt.show()


# In[ ]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_getmm')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("Average metatranscriptome expression (GeTMM)", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
# plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
# plt.ylim(0, 1)

# handles, labels = plt.gca().get_legend_handles_labels()
# plt.legend(reversed(handles), reversed(labels), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_MAG_stacked_area_{treatment_to_plot}.pdf", dpi=300)
plt.show()


# In[ ]:


## Create a separate figure just for the legend
fig_legend = plt.figure(figsize=(8, 6))
# Create proxy artists (patches) for the legend

legend_elements = [Patch(facecolor=mag_color_dict.get(col, 'gray'), 
                        label=col, alpha=0.7) 
                  for col in df_pivot.columns]
# Reverse the order
legend_elements = legend_elements[::-1]
# Create the legend
legend = fig_legend.legend(handles=legend_elements, 
                          loc='center',
                          ncol=1,  # Adjust number of columns as needed
                          fontsize=20)
                          # frameon=True)
                          # fancybox=True,
                          # shadow=True)
# Remove axes
fig_legend.gca().set_axis_off()
# Save the legend
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-B_MAG_legend_only_simplified.pdf",
    dpi=300,
    bbox_inches='tight'
)
plt.show()


# In[ ]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]
# df_treatment = df_treatment.sort_values(by='taxon_grouped_prop_abundance', ascending=False)

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("MAG relative abundance (GeTMM)", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-B_MAG_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[ ]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]
# df_treatment = df_treatment.sort_values(by='taxon_grouped_prop_abundance', ascending=False)

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("MAG relative abundance (GeTMM)", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-B_MAG_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[ ]:


# 1. Pivot the UNFILTERED data to wide format
# (We include 'treatment' and 'day' in the index to preserve them)
df_pivot = by_tax_rank_melted_grouped.pivot(
    index=['treatment', 'day'], 
    columns='taxon_grouped', 
    values='taxon_grouped_prop_abundance'
)

# 2. Reset the index so 'treatment' and 'day' become regular columns
df_pivot = df_pivot.reset_index()

# 3. Create a single-row DataFrame for the colors matching the columns of df_pivot
# We seed it with specific labels for the 'treatment' and 'day' columns
color_row_dict = {'treatment': 'COLOR_CODE', 'day': 'COLOR_CODE'}

# Fill in the rest of the dictionary with the colors for each taxon
for col in df_pivot.columns:
    if col not in ['treatment', 'day']:
        # Get the color from your dict, default to None if missing
        color_row_dict[col] = mag_color_dict.get(col, None)

color_df_row = pd.DataFrame([color_row_dict])

# 4. Concatenate the color row to the TOP of the dataframe
df_pivot = pd.concat([color_df_row, df_pivot], ignore_index=True)
df_pivot.to_csv('data/1B_MAG_prop_geTMM_abundance.csv', index=False)


# In[ ]:


by_tax_rank_melted


# In[ ]:


by_tax_rank_melted.loc[
(by_tax_rank_melted['day'] > 0)
& (by_tax_rank_melted['highest_host_tax_rank'] == 'g__JAEXAI01')
]


# In[ ]:


by_tax_rank_melted.loc[
(by_tax_rank_melted['day'] > 0)
& (by_tax_rank_melted['highest_host_tax_rank'] == 'g__JAATFL01')
]


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 004-MAG-relative-abundance-over-time-getmm.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/01-B_MAG_relative_abundance_over_time_getmm')


# In[ ]:





# In[ ]:




