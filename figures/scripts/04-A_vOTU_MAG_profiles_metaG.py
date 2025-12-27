#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# # The goal is to look at the relative abundance of JAGFXR01 and viral populations within the same samples, and make an estimate on how many viruses to hosts we had

# In[2]:


# host = 'g__Pseudomonas_E'
# host = 'g__Clostridium'
# host = 'g__Paludibacter'
host = "g__JAGFXR01"


# # Host coverage data

# In[3]:


MAG_metaG_meancov = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/coverm/metaG_MAG_mean_97_75_75.tsv', sep='\t')
MAG_metaG_meancov.columns = ['MAG', 'STM_0716_E_M_E001',
       'STM_0716_E_M_E054',
       'STM_0716_E_M_E061',
       'STM_0716_E_M_E034',
       'STM_0716_E_M_E065',
       'STM_0716_E_M_E026',
       'STM_0716_E_M_E058',
       'STM_0716_E_M_E050',
       'STM_0716_E_M_E069',
       'STM_0716_E_M_E030']

MAG_metaG_meancov = MAG_metaG_meancov.melt(id_vars='MAG', var_name='Sample', value_name='meancov')
metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaG_sample_metadata.csv')

MAG_metaG_meancov = MAG_metaG_meancov.merge(metadata, on='Sample', how='left')

# remove CT samples
MAG_metaG_meancov = MAG_metaG_meancov.loc[MAG_metaG_meancov['treatment'] != 'CT']
MAG_metaG_meancov

# filter for active MAGs
MAG_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_abundance_with_metadata.tsv', usecols=['MAG', 'GTDB'], sep='\t')
MAG_metaT = MAG_metaT.drop_duplicates()
MAG_metaG_meancov = MAG_metaT.merge(MAG_metaG_meancov, on='MAG', how='left') # gets rid of MAGs that were not active

MAG_metaG_meancov[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = MAG_metaG_meancov['GTDB'].str.split(';', expand=True)

# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

MAG_metaG_meancov['highest_host_tax_rank'] = MAG_metaG_meancov.apply(get_highest_resolution, axis=1)

MAG_metaG_meancov['day'] = MAG_metaG_meancov['timepoint'].apply(lambda x: int(x.replace('day', '')))
MAG_metaG_meancov = MAG_metaG_meancov.drop(columns=['timepoint'])
MAG_metaG_meancov


# In[4]:


MAG_metaG = MAG_metaG_meancov


# # vOTUs

# In[5]:


vOTU_metaG_meancov = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaG/coverm_trimmed_mean/metaG_vOTU_trimmed_mean_90_75_10.tsv', sep='\t')
vOTU_metaG_meancov.columns = ['vOTU',
                'STM_0716_E_M_E030',
                'STM_0716_E_M_E034',
                'STM_0716_E_M_E069',
                'STM_0716_E_M_E026',
                'STM_0716_E_M_E061',
                'STM_0716_E_M_E058',
                'STM_0716_E_M_E050',
                'STM_0716_E_M_E054',
                'STM_0716_E_M_E065',
                'STM_0716_E_M_E001'
]
vOTU_metaG_meancov


# In[6]:


vOTU_metaG_meancov = vOTU_metaG_meancov.melt(id_vars='vOTU', var_name='Sample', value_name='meancov')
vOTU_metaG_meancov = vOTU_metaG_meancov.merge(metadata, on='Sample', how='left')

# remove CT samples
vOTU_metaG_meancov = vOTU_metaG_meancov.loc[vOTU_metaG_meancov['treatment'] != 'CT']
vOTU_metaG_meancov


# In[7]:


# merge iPHoP predictions
iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/active_vOTUs_active_MAGs_for_cytoscape.tsv', sep='\t')
iphop_df.rename(columns={'Virus': 'vOTU'}, inplace=True)
vOTU_metaG_meancov = vOTU_metaG_meancov.merge(iphop_df[['vOTU', 'highest_host_tax_rank']], on='vOTU', how='left')
vOTU_metaG_meancov


# In[8]:


vOTU_metaG_meancov['day'] = vOTU_metaG_meancov['timepoint'].apply(lambda x: int(x.replace('day', '')))
vOTU_metaG_meancov = vOTU_metaG_meancov.drop(columns=['timepoint'])
vOTU_metaG_meancov


# In[9]:


vOTU_metaG = vOTU_metaG_meancov


# In[10]:


color_dict = {
    '1': '#A6761D',
    '2': '#FFFF99',
    '3': '#A1D99B',
    '4': '#EF3B2C',
    '5': '#1AC7C2',
    
    '2,3,4': '#FCBBA1',
    '2,3,4,5': '#252525',
    '2,4': '#238B45',
    '2,4,5': '#6BAED6',
    '2,5': '#F768A1',
    '3,4': '#807DBA',
    '3,5': '#08589E',

    '': '#808080'
}


# In[11]:


def plot_metaG_profiles(host):

    # filter to host only
    host_df_melted = MAG_metaG.loc[MAG_metaG['highest_host_tax_rank'] == host]
    vOTU_df_melted = vOTU_metaG_meancov.loc[vOTU_metaG_meancov['highest_host_tax_rank'] == host]

    # make color dict from indicator vOTUs
    indicators = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv')
    indicators = indicators.loc[indicators['p_value'] <= 0.01][['vOTU', 'significant_clusters']]
    
    vOTU_df_melted = vOTU_df_melted.merge(indicators, on='vOTU', how='left').fillna('')
    vOTU_df_melted['color'] = vOTU_df_melted['significant_clusters'].apply(lambda x: color_dict[x])
    vOTU_df_melted.rename(columns={'vOTU': 'MAG'}, inplace=True)

    # order shared columns the same
    host_df_melted = host_df_melted[['MAG', 'Sample', 'treatment', 'day', 'meancov']]
    vOTU_df_melted = vOTU_df_melted[['MAG', 'Sample', 'treatment', 'day', 'meancov', 'significant_clusters', 'color']]

    # merge dataframes
    merged_df = pd.concat([host_df_melted, vOTU_df_melted], axis=0)

    # fill NAs
    merged_df['color'] = merged_df['color'].fillna('black')
    merged_df['significant_clusters'] = merged_df['color'].fillna('black')

    # Duplicate day 0 and label catechin
    merged_modified = merged_df.copy()
    
    # Filter for unamended day 0 rows
    unamended_day0 = merged_df[
        (merged_df['treatment'] == 'unamended') & 
        (merged_df['day'] == 0)
    ].copy()
    
    # Change treatment to 'catechin' in the copied rows
    unamended_day0['treatment'] = 'catechin'
    
    # Concatenate the original and modified rows
    merged_df = pd.concat([merged_modified, unamended_day0], ignore_index=True)

    # take mean of each set of sample replicates, should be the same since there are no replicates
    merged_df_mean = merged_df.groupby(['MAG', 'treatment', 'day', 'significant_clusters', 'color']).agg({'meancov': 'mean'}).reset_index()

    # add provirus and MAG designations for labeling
    merged_df_mean['is_provirus'] = merged_df_mean['MAG'].str.contains('provirus')
    merged_df_mean['is_MAG'] = ~merged_df_mean['MAG'].str.contains('contig_')

    # log transform meancov since values are a broad range
    merged_df_mean['log10p_meancov'] = merged_df_mean['meancov'].apply(lambda x: np.log10(x+1))

    # create color dictionary to map colors to correct values
    mag_color_dict = dict(
        zip(
            merged_df_mean.drop_duplicates('MAG')['MAG'],
            merged_df_mean.drop_duplicates('MAG')['color']
        )
    )

    # Set style and figure size
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(
        2, 1, figsize=(6,11), sharex=True
    )  # 2 rows, 1 column
    
    treatments = ['unamended', 'catechin']  # order matters
    
    for ax, treatment in zip(axes, treatments):
        # Filter dataframe by treatment
        treatment_data = merged_df_mean[merged_df_mean['treatment'] == treatment]
        treatment_data = treatment_data.sort_values('is_MAG', ascending=True)
    
        for mag_name in treatment_data['MAG'].unique():
            # Filter data for this specific MAG and sort by day
            mag_data = treatment_data[treatment_data['MAG'] == mag_name].sort_values('day')
    
            # First row for style decisions
            first_row = mag_data.iloc[0]
    
            linestyle = '--' if first_row['is_MAG'] else '-'
            linewidth = 3 if first_row['is_provirus'] or first_row['is_MAG'] else 1
            alpha = 0.7 if first_row['is_MAG'] else 1
    
            # Plot
            ax.plot(
                mag_data['day'],
                mag_data['log10p_meancov'],
                label=mag_name,
                # color=mag_color_dict[mag_name],
                color='black',
                linestyle=linestyle,
                linewidth=linewidth,
                marker='o',
                markersize=4,
                alpha=alpha
            )
    
        # Format axes
        ax.set_title(f"{treatment}", fontsize=12)
        ax.set_xticks([0, 14, 21, 35])
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_ylabel('meancov (log10(x+1))', fontsize=13)
    
    # Add shared x-axis label
    axes[-1].set_xlabel("Day", fontsize=14)
    plt.suptitle(host, fontsize=16, horizontalalignment='right')
    plt.legend().remove()
    plt.tight_layout()
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_{host}_vOTU_metaG.png', dpi=300, bbox_inches='tight')
    plt.show()

    return merged_df


# In[12]:


plot_metaG_profiles('g__Methanothrix')


# In[13]:


plot_metaG_profiles('g__Methanosarcina')


# In[14]:


plot_metaG_profiles('g__Methanobacterium_A')


# In[15]:


plot_metaG_profiles('g__Methanobacterium_B')


# In[16]:


plot_metaG_profiles('g__Paludibacter')


# In[17]:


merged_df = plot_metaG_profiles('g__JAGFXR01')


# In[18]:


merged_df['MAG'].unique()


# In[19]:


plot_metaG_profiles('g__Clostridium')


# In[20]:


plot_metaG_profiles('g__Pseudomonas_E')


# In[21]:


def plot_metaG_profiles_no_color(host):

    # filter to host only
    host_df_melted = MAG_metaG.loc[MAG_metaG['highest_host_tax_rank'] == host]
    vOTU_df_melted = vOTU_metaG_meancov.loc[vOTU_metaG_meancov['highest_host_tax_rank'] == host]

    # make color dict from indicator vOTUs
    indicators = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv')
    indicators = indicators.loc[indicators['p_value'] <= 0.01][['vOTU', 'significant_clusters']]
    
    vOTU_df_melted = vOTU_df_melted.merge(indicators, on='vOTU', how='left').fillna('')
    vOTU_df_melted['color'] = vOTU_df_melted['significant_clusters'].apply(lambda x: color_dict[x])
    vOTU_df_melted.rename(columns={'vOTU': 'MAG'}, inplace=True)

    # order shared columns the same
    host_df_melted = host_df_melted[['MAG', 'Sample', 'treatment', 'day', 'meancov']]
    vOTU_df_melted = vOTU_df_melted[['MAG', 'Sample', 'treatment', 'day', 'meancov', 'significant_clusters', 'color']]

    # merge dataframes
    merged_df = pd.concat([host_df_melted, vOTU_df_melted], axis=0)

    # fill NAs
    merged_df['color'] = merged_df['color'].fillna('#808080')
    merged_df['significant_clusters'] = merged_df['color'].fillna('')

    # Duplicate day 0 and label catechin
    merged_modified = merged_df.copy()
    
    # Filter for unamended day 0 rows
    unamended_day0 = merged_df[
        (merged_df['treatment'] == 'unamended') & 
        (merged_df['day'] == 0)
    ].copy()
    
    # Change treatment to 'catechin' in the copied rows
    unamended_day0['treatment'] = 'catechin'
    
    # Concatenate the original and modified rows
    merged_df = pd.concat([merged_modified, unamended_day0], ignore_index=True)

    # take mean of each set of sample replicates, should be the same since there are no replicates
    merged_df_mean = merged_df.groupby(['MAG', 'treatment', 'day', 'significant_clusters', 'color']).agg({'meancov': 'mean'}).reset_index()

    # add provirus and MAG designations for labeling
    merged_df_mean['is_provirus'] = merged_df_mean['MAG'].str.contains('provirus')
    merged_df_mean['is_MAG'] = ~merged_df_mean['MAG'].str.contains('contig_')

    # log transform meancov since values are a broad range
    merged_df_mean['log10p_meancov'] = merged_df_mean['meancov'].apply(lambda x: np.log10(x+1))

    # create color dictionary to map colors to correct values
    mag_color_dict = dict(
        zip(
            merged_df_mean.drop_duplicates('MAG')['MAG'],
            merged_df_mean.drop_duplicates('MAG')['color']
        )
    )

    # Set style and figure size
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(
        2, 1, figsize=(6,11), sharex=True
    )  # 2 rows, 1 column
    
    treatments = ['unamended', 'catechin']  # order matters
    
    for ax, treatment in zip(axes, treatments):
        # Filter dataframe by treatment
        treatment_data = merged_df_mean[merged_df_mean['treatment'] == treatment]
    
        for mag_name in treatment_data['MAG'].unique():
            # Filter data for this specific MAG and sort by day
            mag_data = treatment_data[treatment_data['MAG'] == mag_name].sort_values('day')
    
            # First row for style decisions
            first_row = mag_data.iloc[0]
    
            linestyle = '--' if first_row['is_MAG'] else '-'
            linewidth = 3 if first_row['is_provirus'] or first_row['is_MAG'] else 2
            alpha = 1 if first_row['is_MAG'] else 0.5
    
            # Plot
            ax.plot(
                mag_data['day'],
                mag_data['log10p_meancov'],
                label=mag_name,
                color="black",
                linestyle=linestyle,
                linewidth=linewidth,
                marker='o',
                markersize=4,
                alpha=alpha
            )
    
        # Format axes
        ax.set_title(f"{treatment}", fontsize=12)
        ax.set_xticks([0, 14, 21, 35])
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_ylabel('meancov (log10(x+1))', fontsize=13)
    
    # Add shared x-axis label
    axes[-1].set_xlabel("Day", fontsize=14)
    plt.suptitle(host, fontsize=16, horizontalalignment='right')
    plt.legend().remove()
    plt.tight_layout()
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_{host}_vOTU_metaG_no_color.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    return merged_df


# In[22]:


plot_metaG_profiles_no_color('g__JAGFXR01')


# In[23]:


get_ipython().system('jupyter nbconvert --to script 007-vOTU_vs_host_metaG.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/03-A_vOTU_MAG_profiles_metaG')


# In[ ]:




