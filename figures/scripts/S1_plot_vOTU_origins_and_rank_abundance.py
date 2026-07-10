#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

import matplotlib as mpl
from print_versions import print_versions

# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'

print_versions(globals())


# In[2]:


filtered_vOTUs_metadata_file = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv'

with open(filtered_vOTUs_metadata_file, 'r') as f:
    filtered_vOTUs = [i.split('\t')[0] for i in f.readlines()]
    filtered_vOTUs = list(set(filtered_vOTUs[1:]))


# In[3]:


len(filtered_vOTUs)


# In[4]:


vOTUs_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.tsv', header=None, sep='\t')
vOTUs_df.columns = ['rep', 'members']
filtered_vOTUs_df = vOTUs_df.loc[vOTUs_df['rep'].isin(filtered_vOTUs)].copy()
filtered_vOTUs_df


# In[5]:


# get header IDs from /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/grantham_assemblies_5kb
SUN_START=1166228
SUN_END=1171278
sun2024_contigs_list = ['contig_' + str(i) for i in list(range(SUN_START, SUN_END+1))]

EMERSON_START=1171279
EMERSON_END=1173181
emerson2018_contigs_list = ['contig_' + str(i) for i in list(range(EMERSON_START, EMERSON_END+1))]

GRANTHAM_START=1
GRANTHAM_END=1166227
grantham_contigs_list = ['contig_' + str(i) for i in list(range(GRANTHAM_START, GRANTHAM_END+1))]


# In[6]:


def is_in_cluster(lst, df):
    """
    Marks '1' in the 'genomad' column if at least one member of 'lst' 
    is found in the 'members' column of df, otherwise marks '0'.

    Parameters:
        lst (list): List of contig members to check.
        df (pd.DataFrame): DataFrame with 'members' column (comma-separated values).

    Returns:
        pd.Series: Updated 'genomad' column with 1s and 0s.

    Example:
        lst = ['a', 'd']
        df:
            rep  members  genomad
            a    a        0
            b    b,d      0
            c    c,e,f    0
        
        df['genomad'] = is_in_cluster(lst, df)
        
        Output:
            rep  members  genomad
            a    a        1
            b    b,d      1
            c    c,e,f    0  
    """
    vOTU_set = set(lst)
    
    def check_members(cluster):
        # Clean each member by removing everything after '|'
        members = [m.split('|')[0] for m in cluster.split(',')]
        return int(bool(set(members) & vOTU_set))

    return df['members'].apply(check_members)


# In[7]:


filtered_vOTUs_df.loc[:, 'grantham'] = is_in_cluster(grantham_contigs_list, filtered_vOTUs_df)
filtered_vOTUs_df.loc[:, 'sun2024'] = is_in_cluster(sun2024_contigs_list, filtered_vOTUs_df)
filtered_vOTUs_df.loc[:, 'emerson2018'] = is_in_cluster(emerson2018_contigs_list, filtered_vOTUs_df)
filtered_vOTUs_df.head()


# In[8]:


filtered_vOTUs_df.to_csv('data/S1A_votu_reference_database_origins.csv')


# In[9]:


grantham_set = set(filtered_vOTUs_df.loc[filtered_vOTUs_df['grantham'] >= 1]['rep'])
sun2024_set = set(filtered_vOTUs_df.loc[filtered_vOTUs_df['sun2024'] >= 1]['rep'])
emerson2018_set = set(filtered_vOTUs_df.loc[filtered_vOTUs_df['emerson2018'] >= 1]['rep'])


# In[10]:


len(grantham_set)


# In[11]:


len(sun2024_set)


# In[12]:


len(emerson2018_set)


# In[13]:


# Define colorblind-friendly colors (Paul Tol's scheme)
# colors = ["#0072B2", "#E69F00", "#009E73"]  # Blue, Orange, Green
colors = ["white", "white", "white"]
# Create a high-resolution figure
plt.figure(figsize=(8,6), dpi=300)

# Create Venn diagram with custom colors
venn = venn3([grantham_set, sun2024_set, emerson2018_set], 
             set_labels=('This study', 'Sun & Pratama et al. 2024', 'Emerson et al. 2018'),
             set_colors=colors, alpha=0.7)  # Adjust transparency for better contrast

# Customize borders
for patch in venn.patches:
    if patch:
        patch.set_edgecolor("black")  # Black borders
        patch.set_linewidth(1.5)  # Thicker borders

# Improve text visibility
for label in venn.subset_labels:
    if label:
        label.set_fontsize(12)  # Increase font size
        # label.set_fontweight('bold')  # Make text bold
        label.set_bbox(dict(facecolor='lightgrey', edgecolor='none', boxstyle='round,pad=0.2'))  # White background

# save the plot
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1A_vOTU-origins-no-cutoff.pdf', dpi=300)

# Show the plot
plt.show()


# In[14]:


9562 / (9562-3534)


# ### Filter for >10kb to compare

# In[15]:


vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', usecols=['vOTU', 'length'], sep='\t').drop_duplicates(subset='vOTU')
vOTU_metadata.head()


# In[16]:


filtered_vOTUs_10kb_df = filtered_vOTUs_df.merge(vOTU_metadata, left_on='rep', right_on='vOTU', how='left')
filtered_vOTUs_10kb_df = filtered_vOTUs_10kb_df.loc[filtered_vOTUs_10kb_df['length'] >= 10000]
filtered_vOTUs_10kb_df = filtered_vOTUs_10kb_df.drop(columns=['vOTU', 'length'])
filtered_vOTUs_10kb_df


# In[17]:


grantham_10kb_set = set(filtered_vOTUs_10kb_df.loc[(filtered_vOTUs_10kb_df['grantham'] >= 1)]['rep'])
sun2024_10kb_set = set(filtered_vOTUs_10kb_df.loc[filtered_vOTUs_10kb_df['sun2024'] >= 1]['rep'])
emerson2018_10kb_set = set(filtered_vOTUs_10kb_df.loc[filtered_vOTUs_10kb_df['emerson2018'] >= 1]['rep'])


# In[18]:


# Define colorblind-friendly colors (Paul Tol's scheme)
# colors = ["#0072B2", "#E69F00", "#009E73"]  # Blue, Orange, Green

# Create a high-resolution figure
plt.figure(figsize=(8,6), dpi=300)

# Create Venn diagram with custom colors
venn = venn3([grantham_10kb_set, sun2024_10kb_set, emerson2018_10kb_set], 
             set_labels=('This study', 'Sun & Pratama et al. 2024', 'Emerson et al. 2018'),
             set_colors=colors, alpha=0.7)  # Adjust transparency for better contrast

# Customize borders
for patch in venn.patches:
    if patch:
        patch.set_edgecolor("black")  # Black borders
        patch.set_linewidth(1.5)  # Thicker borders

# Improve text visibility
for label in venn.subset_labels:
    if label:
        label.set_fontsize(12)  # Increase font size
        # label.set_fontweight('bold')  # Make text bold
        label.set_bbox(dict(facecolor='lightgrey', edgecolor='none', boxstyle='round,pad=0.2'))  # White background

# Adjust title
plt.title("Venn Diagram of Stordalen vOTU Reference Database", fontsize=14, fontweight='bold')

# save the plot
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX-vOTU-origins-no-cutoff-10kb.pdf', dpi=300)

# Show the plot
plt.show()


# # The below is run after identifying which vOTUs are active (001-identify-active-vOTU-and-calc-getmm-abundance.ipynb)

# In[19]:


# filter for only active vOTUs
getmms_vOTU_corrected_wide = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
active_vOTUs = list(getmms_vOTU_corrected_wide.loc[getmms_vOTU_corrected_wide.iloc[:,1:27].sum(axis=1) > 0].vOTU.unique())
active_vOTUs[:5]


# In[20]:


active_filtered_vOTUs_df = filtered_vOTUs_df.loc[filtered_vOTUs_df['rep'].isin(active_vOTUs)]

active_grantham_set = set(active_filtered_vOTUs_df.loc[active_filtered_vOTUs_df['grantham'] >= 1]['rep'])
active_sun2024_set = set(active_filtered_vOTUs_df.loc[active_filtered_vOTUs_df['sun2024'] >= 1]['rep'])
active_emerson2018_set = set(active_filtered_vOTUs_df.loc[active_filtered_vOTUs_df['emerson2018'] >= 1]['rep'])


# In[21]:


active_filtered_vOTUs_df.to_csv('data/S1A_votu_active_origins.csv')


# In[22]:


len(active_grantham_set)


# In[23]:


# Define colorblind-friendly colors (Paul Tol's scheme)
# colors = ["#0072B2", "#E69F00", "#009E73"]  # Blue, Orange, Green

# Create a high-resolution figure
plt.figure(figsize=(8,6), dpi=300)

# Create Venn diagram with custom colors
venn = venn3([active_grantham_set, active_sun2024_set, active_emerson2018_set], 
             set_labels=('This study', 'Sun & Pratama 2024', 'Emerson 2018'),
             set_colors=colors, alpha=0.7)  # Adjust transparency for better contrast

# Customize borders
for patch in venn.patches:
    if patch:
        patch.set_edgecolor("black")  # Black borders
        patch.set_linewidth(1.5)  # Thicker borders

# Improve text visibility
for label in venn.subset_labels:
    if label:
        label.set_fontsize(12)  # Increase font size
        # label.set_fontweight('bold')  # Make text bold
        label.set_bbox(dict(facecolor='lightgrey', edgecolor='none', boxstyle='round,pad=0.2'))  # White background

# Adjust set label positions manually to prevent overlap
if venn.set_labels:
    venn.set_labels[0].set_position((-0.5, 0.4))  # Move "Grantham"
    venn.set_labels[1].set_position((0.45, 0.25))   # Move "Sun 2024"
    venn.set_labels[2].set_position((0.6, -0.30))    # Move "Emerson 2018" downward slightly

# Adjust title
# plt.title("Venn Diagram of active vOTU Overlaps", fontsize=14, fontweight='bold')

plt.tight_layout()
# save the plot
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1A_active_vOTU_origins.pdf', dpi=300)

# Show the plot
plt.show()


# # Check rank abundance to make sure it's not a sequencing depth artifact, and that the new vOTUs are not just rares we are getting from deeper sequencing

# In[24]:


abundance_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
abundance_df


# In[25]:


abundance_long = abundance_df.melt(id_vars='vOTU', var_name='Sample', value_name='abundance')
abundance_long.head()


# In[26]:


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
    'STM_0716_E_M_E121_unamended_07': 1,
       'STM_0716_E_M_E122_unamended_07': 2,
    'STM_0716_E_M_E123_unamended_07': 3,
       'STM_0716_E_M_E125_CT_07': 1,
    'STM_0716_E_M_E126_CT_07': 2,
       'STM_0716_E_M_E127_CT_07': 3,
    'STM_0716_E_M_E129_catechin_07': 1,
       'STM_0716_E_M_E130_catechin_07': 2,
    'STM_0716_E_M_E131_catechin_07': 3
}.items())

replicate_frame['Sample'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[0])
replicate_frame['treatment'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[1])
replicate_frame['day'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[2])
replicate_frame.columns = ['sample_treatment_day', 'replicate', 'Sample', 'treatment', 'day']
replicate_frame['treatment_day_replicate'] = replicate_frame['treatment'] + '_' + replicate_frame['day'].astype(str) + '_' + replicate_frame['replicate'].astype(str)
replicate_frame.head()


# In[27]:


abundance_long = abundance_long.merge(replicate_frame[['Sample', 'treatment']], on ='Sample', how='left')


# In[28]:


filtered_vOTUs_df.columns = ['vOTU', 'members', 'grantham', 'sun2024', 'emerson2018']


# In[29]:


abundance_long = abundance_long.merge(filtered_vOTUs_df, on='vOTU', how='left')
abundance_long = abundance_long.dropna(subset='grantham')


# In[30]:


abundance_agg = abundance_long.groupby(['vOTU', 'grantham', 'sun2024', 'emerson2018', 'treatment']).agg({'abundance': 'max'})


# In[31]:


abundance_agg = abundance_agg.reset_index()


# In[32]:


abundance_agg


# In[33]:


df_sorted = abundance_agg.sort_values(by='abundance', ascending=False).reset_index(drop=True)
df_sorted['rank'] = df_sorted.index + 1


# In[34]:


df_sorted['log10_abundance'] = df_sorted['abundance'].apply(lambda x: np.log10(x + 1))


# In[35]:


# make sure all vOTUs have an origin assignment.
df_sorted.loc[df_sorted[['grantham', 'sun2024', 'emerson2018']].sum(axis=1) == 0]


# In[36]:


df_sorted


# In[37]:


df_sorted.loc[(df_sorted['grantham'] + df_sorted['sun2024'] + df_sorted['emerson2018']) == 0]


# # Get rank for only active

# In[38]:


vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
provirus_df = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()


# In[39]:


df_sorted = df_sorted.merge(provirus_df, on='vOTU', how='left')


# In[40]:


membership_bool = df_sorted[['grantham', 'sun2024', 'emerson2018']].T.astype(bool)


# ## Sort by provirus or not

# In[41]:


df_sorted = df_sorted.sort_values(by=['is_provirus_aggregated', 'log10_abundance'], ascending=False)
# Update rank after sorting
df_sorted['rank'] = np.arange(len(df_sorted))


# In[42]:


df_sorted.to_csv('S1D_active_votu_rank_abundance_source_data.csv


# In[40]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Your existing aggregation code
abundance_agg = abundance_long.groupby(['vOTU', 'grantham', 'sun2024', 'emerson2018', 'treatment']).agg({'abundance': 'max'})
abundance_agg = abundance_agg.reset_index()

# Pivot to get separate columns for each treatment
abundance_pivot = abundance_agg.pivot_table(
    index=['vOTU', 'grantham', 'sun2024', 'emerson2018'], 
    columns='treatment', 
    values='abundance', 
    fill_value=0
).reset_index()

# Rank based on catechin treatment abundance
abundance_pivot = abundance_pivot.sort_values(by='unamended', ascending=False).sort_values(by='catechin', ascending=False).reset_index(drop=True)
abundance_pivot['rank'] = abundance_pivot.index + 1

# Calculate log10 abundance for both treatments
abundance_pivot['log10_catechin'] = abundance_pivot['catechin'].apply(lambda x: np.log10(x + 1))
abundance_pivot['log10_unamended'] = abundance_pivot['unamended'].apply(lambda x: np.log10(x + 1))

# Merge with provirus metadata
vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
provirus_df = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()
abundance_pivot = abundance_pivot.merge(provirus_df, on='vOTU', how='left')

# Sort by provirus status first, then by catechin abundance
abundance_pivot = abundance_pivot.sort_values(by=['is_provirus_aggregated', 'catechin'], ascending=False)
abundance_pivot['rank'] = np.arange(len(abundance_pivot))

# Create membership boolean matrix for heatmap
membership_bool = abundance_pivot[['grantham', 'sun2024', 'emerson2018']].T.astype(bool)

# Create RGB color matrix for heatmap
color_matrix = np.zeros((membership_bool.shape[0], membership_bool.shape[1], 3))  # RGB
for i in range(membership_bool.shape[0]):
    for j in range(membership_bool.shape[1]):
        if membership_bool.iloc[i, j]:  # origin value == 1
            if abundance_pivot.iloc[j]['is_provirus_aggregated']:  # is_provirus == True
                color_matrix[i, j] = [1, 0, 0]  # Red
            else:
                color_matrix[i, j] = [0, 0, 0]  # Black
        else:
            color_matrix[i, j] = [1, 1, 1]  # White

# Find the boundary between provirus and non-provirus for vertical line
provirus_boundary = abundance_pivot['is_provirus_aggregated'].sum() - 0.5

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 6),
                               gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

# Rank-abundance plot with two points per vOTU
ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_catechin'], 
            c='#FC9B2D', 
            marker='o', s=20, alpha=0.7, label='Catechin')

ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_unamended'], 
            c='#7ACBC3', 
            marker='o', s=20, alpha=0.7, label='Unamended')

# Add vertical line to separate provirus from non-provirus
ax1.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8, 
            label='Provirus/Non-provirus boundary')

# Add "provirus" text to the top-left corner, left of the vertical line
ax1.text(provirus_boundary * 0.05, 
         ax1.get_ylim()[1] * 0.95, 
         'provirus', 
         ha='left', va='top', fontsize=14, fontweight='bold', color='red')

ax1.set_ylabel('log10(GeTMM + 1)', size=14)
ax1.legend()


# Heatmap with interpolation disabled
ax2.imshow(color_matrix, aspect='auto', interpolation='none')
ax2.set_ylabel('Source DB', size=12)
ax2.set_yticks(range(len(membership_bool.index)))
ax2.set_yticklabels(['This study', 'Sun2024', 'Emerson2018'], size=12)
ax2.set_xticks([])
ax2.set_xlabel('vOTUs (ranked by maximum GeTMM abundance in a catechin sample)', size=14)

# Add vertical line to heatmap as well
ax2.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1D_rank_abundance_treatment_comparison.pdf', dpi=300)
plt.show()


# In[41]:


get_ipython().system('jupyter nbconvert --to script 000-compare_vOTU_origins.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S1_plot_vOTU_origins_and_rank_abundance')


# In[ ]:





# In[ ]:




