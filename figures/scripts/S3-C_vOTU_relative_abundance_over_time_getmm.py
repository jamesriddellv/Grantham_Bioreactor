#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[3]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_vOTUs_wide.tsv', sep='\t')
df = df.set_index('vOTU')
df.head()


# In[5]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.sort_values(by=['vOTU', 'Confidence score']).drop_duplicates('vOTU')
# due to changes in taxonomy, we will change g__Lachnoclostridium to g__Clostridium
# iphop_df['G'] = iphop_df['G'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
# iphop_df['highest_host_tax_rank'] = iphop_df['highest_host_tax_rank'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
iphop_df.head()


# In[6]:


# here we correctly handle situations where vOTUs might have had an iphop assignment, but wasn't an active MAG, therefore is called "no_assignment"
df = df.merge(iphop_df, on='vOTU', how='left')


# In[7]:


df


# In[17]:


df['highest_host_tax_rank'] = df['highest_host_tax_rank'].fillna('no_assignment')
sum(df['highest_host_tax_rank'] != 'no_assignment') / len(df)


# In[22]:


keep_cols = ['highest_host_tax_rank', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']


# In[23]:


# Sum getmm values of vOTUs in the same highest host tax rank
by_tax_rank = df.groupby('highest_host_tax_rank').sum().reset_index()[keep_cols]


# In[24]:


by_tax_rank_melted = by_tax_rank.melt(id_vars=['highest_host_tax_rank'], var_name='Sample', value_name='abundance')
by_tax_rank_melted


# In[25]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# In[26]:


by_tax_rank_melted = by_tax_rank_melted.merge(replicate_frame, on='Sample', how='left')


# In[27]:


by_tax_rank_melted


# In[28]:


# Take the mean across replicates, since there are fewer replicates of unamended day 14
by_tax_rank_melted = by_tax_rank_melted.groupby(['highest_host_tax_rank', 'treatment', 'day']).agg({'abundance': 'mean'}).reset_index()


# In[29]:


by_tax_rank_melted


# In[30]:


# Add prop abundance
by_sample_sum_abundance = by_tax_rank_melted.groupby(['treatment', 'day']).agg({'abundance': 'sum'}).reset_index().rename(columns={'abundance': 'total_getmm'})
by_sample_sum_abundance


# In[31]:


by_tax_rank_melted = by_tax_rank_melted.merge(by_sample_sum_abundance, on=['treatment', 'day'], how='left')


# In[32]:


by_tax_rank_melted['prop_abundance'] = by_tax_rank_melted['abundance'] / by_tax_rank_melted['total_getmm']


# In[33]:


by_tax_rank_melted


# In[34]:


by_tax_rank_melted['day'] = by_tax_rank_melted['day'].astype(int)


# In[35]:


# Step 1: Make sure 'taxon' is the correct column name

# Step 2: Group by treatment_day_replicate and taxon, calculate mean prop_abundance
# grouped = (
#     by_tax_rank_melted.groupby(['treatment_day_replicate', 'highest_host_tax_rank_appended'])['prop_abundance']
#     .mean()
#     .reset_index()
# )

# Step 3: Find max prop_abundance per taxon across all groups
max_abundance = by_tax_rank_melted.groupby('highest_host_tax_rank')['prop_abundance'].max()

# Step 4: Identify taxa to keep (those with at least one group ≥ 5%)
taxa_to_keep = max_abundance[max_abundance >= 0.04].index

print(taxa_to_keep)


# In[37]:


taxa_to_keep = ['Other',
 'g__Clostridium',
 'g__Clostridium__contig_507316',
 'g__Clostridium__contig_588865',
 'g__Clostridium__contig_755337',
 'g__Fen-1039',
 'g__Fen-1308',
 'g__JAGFXR01',
 'g__JAGFXR01__contig_591846',
 'g__PALSA-129',
 'g__Paludibacter',
 'g__Paludibacter__contig_724521',
 'g__Pseudomonas_E',
 'g__Terracidiphilus',
 'g__UBA1794',
 'no_assignment',
 'no_assignment__contig_1022629',
 'no_assignment__contig_1093732',
 'no_assignment__contig_709158',
 'no_assignment__contig_713418',
 'no_assignment__contig_773239',
 'no_assignment__contig_776767',
 'no_assignment__contig_861080']

# Step 5: Create a new column in original df with taxon or "other"
by_tax_rank_melted['taxon_grouped'] = by_tax_rank_melted['highest_host_tax_rank'].where(by_tax_rank_melted['highest_host_tax_rank'].isin(taxa_to_keep), 'Other')


# In[38]:


by_tax_rank_melted.columns


# In[39]:


by_tax_rank_melted


# In[40]:


by_tax_rank_melted_grouped = by_tax_rank_melted.groupby(['treatment', 'day', 'taxon_grouped']).agg(taxon_grouped_abundance=('abundance', 'sum'), taxon_grouped_prop_abundance=('prop_abundance', 'sum')).reset_index()
by_tax_rank_melted_grouped


# In[41]:


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


# In[42]:


by_tax_rank_melted_grouped


# In[43]:


mag_color_dict = {
    'Other': 'black',
    
    # g__Clostridium group (red shades) - maintaining existing colors
    'g__Clostridium': '#9E0B0F',                    # Original sangria red
    'g__Clostridium__contig_507316': '#DC143C',     # Crimson (from original)
    'g__Clostridium__contig_588865': '#FF4500',     # Orange red (new - similar shade)
    'g__Clostridium__contig_755337': '#B22222',     # Fire brick red (from original)
    
    # g__Fen-1039 group (orange shades)
    'g__Fen-1039': '#FF8C00',                       # Original dark orange
    
    # g__JAGFXR01 group (purple shades) - maintaining existing colors
    'g__JAGFXR01': '#845EC2',                       # Original purple
    'g__JAGFXR01__contig_591846': '#5C4188',        # Darker purple (from original)
    
    # g__PALSA-129 group (blue shades) - new from second list
    'g__PALSA-129': '#004E89',                      # Original deep blue
    
    # g__Methanosarcina - from first list
    'g__Methanosarcina': '#228B22',                 # Forest green (from original)
    
    # g__Paludibacter group (bright red shades) - maintaining existing colors
    'g__Paludibacter': '#FF1744',                   # Original bright red
    'g__Paludibacter__contig_724521': '#E91E63',    # Pink-red (from original)
    
    # g__Pseudomonas_E group (gold shades)
    'g__Pseudomonas_E': '#FFD700',                  # Original gold
    
    # g__Terracidiphilus group (green shades)
    'g__Terracidiphilus': '#32CD32',                # Original lime green
    
    # g__UBA1794 group (pink shades)
    'g__UBA1794': '#FF69B4',                        # Original hot pink
    
    # no_assignment group (gray shades) - maintaining existing colors
    'no_assignment__contig_1022629': '#C0C0C0',     # Silver (from original)
    'no_assignment__contig_1093732': '#D3D3D3',     # Light gray (from original)
    'no_assignment__contig_709158': '#696969',      # Dim gray (from original)
    'no_assignment__contig_713418': '#778899',      # Light slate gray (from original)
    'no_assignment__contig_773239': '#708090',      # Slate gray (from original)
    'no_assignment__contig_776767': '#556B2F',      # Dark olive green (from original)
    'no_assignment__contig_861080': '#2F4F4F',      # Dark slate gray (from original)
    'no_assignment': '#808080',                     # Medium gray
}


# In[44]:


list(mag_color_dict.keys())


# In[45]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_abundance')
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
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

# Set yticks with custom formatting
import matplotlib.ticker as ticker

# Define the tick values and labels
yticks_values = [0, 5000, 10000, 15000, 20000]
yticks_labels = ['0', '5k', '10k', '15k', '20k']

# Set the yticks
plt.yticks(yticks_values, yticks_labels, fontsize=12)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3-C_vOTU_stacked_area_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[46]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_abundance')
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
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

# Set yticks with custom formatting
import matplotlib.ticker as ticker

# Define the tick values and labels
yticks_values = [0, 5000, 10000, 15000, 20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000]
yticks_labels = ['0', '5k', '10k', '15k', '20k', '40k', '60k', '80k', '100k', '120k', '140k', '160k']

# Set the yticks
plt.yticks(yticks_values, yticks_labels, fontsize=12)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3-C_vOTU_stacked_area_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[50]:


## Create a separate figure just for the legend
fig_legend = plt.figure(figsize=(8, 6))
# Create proxy artists (patches) for the legend
from matplotlib.patches import Patch
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
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_legend_only_simplified.pdf",
    dpi=300,
    bbox_inches='tight'
)
plt.show()


# In[51]:


df.loc[df['highest_host_tax_rank'] == 'g__JAGFXR01']


# In[52]:


102652.241407 / 750.845462


# In[53]:


by_tax_rank_melted_grouped


# In[54]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

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
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[55]:


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
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# # What proportion of JAGFXR01 is contig_591846 activity?

# In[57]:


jag_df = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == 'catechin']


# In[59]:


jag_df = jag_df.loc[jag_df['taxon_grouped'].str.contains('JAGFXR01')]
jag_df


# In[63]:


jag_df = jag_df.merge(jag_df.groupby(['treatment', 'day']).agg({'taxon_grouped_abundance': 'sum'}).reset_index(), on=['treatment','day'], suffixes=['', '_total'])
             


# In[66]:


jag_df['prop_total'] = jag_df['taxon_grouped_abundance'] / jag_df['taxon_grouped_abundance_total']
jag_df


# In[69]:


jag_df.loc[jag_df['taxon_grouped'] == 'g__JAGFXR01__contig_591846']['prop_total'].mean()


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 004-vOTU-relative-abundance-over-time-getmm.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/02-C_vOTU_relative_abundance_over_time_getmm')

