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
host = 'g__Paludibacter'


# # Load metaT datasets

# In[3]:


MAG_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_GTDB_getmm.tsv', sep='\t')
MAG_metaT.head()


# In[4]:


vOTU_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/vOTU_metadata_getmm.tsv', sep='\t')
vOTU_metaT.head()


# In[5]:


# indicator colors from cytoscape
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


# In[6]:


def plot_metaT_profiles(host):

    # filter to host only
    host_df = MAG_metaT.loc[MAG_metaT['highest_host_tax_rank'] == host]
    host_df_melted = host_df.iloc[:,:27].melt(id_vars='MAG', var_name='Sample', value_name='getmm')

    vOTU_df = vOTU_metaT.loc[vOTU_metaT['G'] == host]
    vOTU_df_melted = vOTU_df.iloc[:,:27].melt(id_vars='vOTU', var_name='Sample', value_name='getmm')
    
    # load metadata
    metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv',
                      usecols=['Sample', 'treatment', 'timepoint'])
    metadata['day'] = metadata['timepoint'].apply(lambda x: int(x.replace('day', '')))
    metadata = metadata.drop(columns=['timepoint'])

    # merge to host and vOTU dfs
    host_df_melted = host_df_melted.merge(metadata, on='Sample', how='left')
    vOTU_df_melted = vOTU_df_melted.merge(metadata, on='Sample', how='left')

    # make color dict from indicator vOTUs
    indicators = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv')
    indicators = indicators.loc[indicators['p_value'] <= 0.01][['vOTU', 'significant_clusters']]
    
    vOTU_df_melted = vOTU_df_melted.merge(indicators, on='vOTU', how='left').fillna('')
    vOTU_df_melted.columns = ['MAG', 'Sample', 'getmm', 'treatment', 'day', 'significant_clusters']
    vOTU_df_melted['color'] = vOTU_df_melted['significant_clusters'].apply(lambda x: color_dict[x])

    # order shared columns the same
    host_df_melted = host_df_melted[['MAG', 'Sample', 'treatment', 'day', 'getmm']]
    vOTU_df_melted = vOTU_df_melted[['MAG', 'Sample', 'treatment', 'day', 'getmm', 'significant_clusters', 'color']]

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

    # take mean of each set of sample replicates
    merged_df_mean = merged_df.groupby(['MAG', 'treatment', 'day', 'significant_clusters', 'color']).agg({'getmm': 'mean'}).reset_index()

    # add provirus and MAG designations for labeling
    merged_df_mean['is_provirus'] = merged_df_mean['MAG'].str.contains('provirus')
    merged_df_mean['is_MAG'] = ~merged_df_mean['MAG'].str.contains('contig_')

    # log transform GeTMM since values are a broad range
    merged_df_mean['log10p_getmm'] = merged_df_mean['getmm'].apply(lambda x: np.log10(x+1))

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
        # Make a twin axis for MAGs
        ax_mag = ax.twinx()
    
        # Filter dataframe by treatment
        treatment_data = merged_df_mean[merged_df_mean['treatment'] == treatment]
    
        for mag_name in treatment_data['MAG'].unique():
            # Filter data for this specific MAG and sort by day
            mag_data = treatment_data[treatment_data['MAG'] == mag_name].sort_values('day')
    
            # First row for style decisions
            first_row = mag_data.iloc[0]
    
            linestyle = '--' if first_row['is_MAG'] else '-'
            linewidth = 3 if first_row['is_provirus'] or first_row['is_MAG'] else 1
            alpha = 1 if first_row['significant_clusters'] != '' else 0.2
    
            # Choose which axis to plot on
            target_ax = ax_mag if first_row['is_MAG'] else ax
    
            target_ax.plot(
                mag_data['day'],
                mag_data['log10p_getmm'],
                label=mag_name,
                color=mag_color_dict[mag_name],
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
        ax.set_ylabel('vOTU GeTMM (log10(x+1))', fontsize=13)
    
        # Label MAG axis separately
        ax_mag.set_ylabel("MAG GeTMM (log10(x+1))", fontsize=13, color="gray")
        ax_mag.tick_params(axis='y', colors="gray")
    
    # Shared x-axis label
    axes[-1].set_xlabel("Day", fontsize=14)
    
    plt.suptitle(host, fontsize=16, ha="right")
    plt.tight_layout()
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_{host}_vOTU_metaT.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    return merged_df


# In[7]:


def plot_metaT_profiles_no_MAG(host):
    # filter to host only
    host_df = MAG_metaT.loc[MAG_metaT['highest_host_tax_rank'] == host]
    host_df_melted = host_df.iloc[:,:27].melt(id_vars='MAG', var_name='Sample', value_name='getmm')

    vOTU_df = vOTU_metaT.loc[vOTU_metaT['G'] == host]
    vOTU_df_melted = vOTU_df.iloc[:,:27].melt(id_vars='vOTU', var_name='Sample', value_name='getmm')
    
    # load metadata
    metadata = pd.read_csv(
        '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv',
        usecols=['Sample', 'treatment', 'timepoint']
    )
    metadata['day'] = metadata['timepoint'].apply(lambda x: int(x.replace('day', '')))
    metadata = metadata.drop(columns=['timepoint'])

    # merge metadata
    host_df_melted = host_df_melted.merge(metadata, on='Sample', how='left')
    vOTU_df_melted = vOTU_df_melted.merge(metadata, on='Sample', how='left')

    # indicator colors
    indicators = pd.read_csv(
        '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv'
    )
    indicators = indicators.loc[indicators['p_value'] <= 0.01][['vOTU', 'significant_clusters']]
    
    vOTU_df_melted = vOTU_df_melted.merge(indicators, on='vOTU', how='left').fillna('')
    vOTU_df_melted.columns = ['MAG', 'Sample', 'getmm', 'treatment', 'day', 'significant_clusters']
    vOTU_df_melted['color'] = vOTU_df_melted['significant_clusters'].apply(lambda x: color_dict[x])

    # combine
    merged_df = pd.concat([host_df_melted, vOTU_df_melted], axis=0)
    merged_df['color'] = merged_df['color'].fillna('#808080')
    merged_df['significant_clusters'] = merged_df['color'].fillna('')

    # duplicate unamended day 0 as catechin day 0
    merged_modified = merged_df.copy()
    unamended_day0 = merged_df[
        (merged_df['treatment'] == 'unamended') & 
        (merged_df['day'] == 0)
    ].copy()
    unamended_day0['treatment'] = 'catechin'
    merged_df = pd.concat([merged_modified, unamended_day0], ignore_index=True)

    # aggregate replicates
    merged_df_mean = merged_df.groupby(
        ['MAG', 'treatment', 'day', 'significant_clusters', 'color']
    ).agg({'getmm': 'mean'}).reset_index()

    # annotate and transform
    merged_df_mean['is_provirus'] = merged_df_mean['MAG'].str.contains('provirus')
    merged_df_mean['is_MAG'] = ~merged_df_mean['MAG'].str.contains('contig_')
    merged_df_mean['log10p_getmm'] = merged_df_mean['getmm'].apply(lambda x: np.log10(x + 1))

    # filter: only non-MAGs
    merged_df_mean = merged_df_mean[merged_df_mean['is_MAG'] == False]

    # color mapping
    mag_color_dict = dict(
        zip(
            merged_df_mean.drop_duplicates('MAG')['MAG'],
            merged_df_mean.drop_duplicates('MAG')['color']
        )
    )

    # plotting
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(2, 1, figsize=(6, 11), sharex=True)

    treatments = ['unamended', 'catechin']

    for ax, treatment in zip(axes, treatments):
        treatment_data = merged_df_mean[merged_df_mean['treatment'] == treatment]
        
        for mag_name in treatment_data['MAG'].unique():
            mag_data = treatment_data[treatment_data['MAG'] == mag_name].sort_values('day')

            first_row = mag_data.iloc[0]
            linewidth = 3 if first_row['is_provirus'] else 1
            alpha = 1 if first_row['significant_clusters'] != '' else 0.2

            ax.plot(
                mag_data['day'],
                mag_data['log10p_getmm'],
                label=mag_name,
                # color=mag_color_dict[mag_name],
                color='black',
                linestyle='-',
                linewidth=linewidth,
                marker='o',
                markersize=4,
                alpha=alpha
            )

        ax.set_title(f"{treatment}", fontsize=12)
        ax.set_xticks([0, 14, 21, 35])
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_ylabel('vOTU GeTMM (log10(x+1))', fontsize=13)

    axes[-1].set_xlabel("Day", fontsize=14)
    plt.suptitle(host, fontsize=16, ha="right")
    plt.tight_layout()
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_{host}_vOTU_metaT_no_MAG.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    return merged_df_mean


# In[8]:


jag = plot_metaT_profiles_no_MAG('g__JAGFXR01')


# In[9]:


plot_metaT_profiles('g__Paludibacter')


# In[10]:


jag = plot_metaT_profiles('g__JAGFXR01')


# In[11]:


plot_metaT_profiles('g__Clostridium')


# In[12]:


plot_metaT_profiles('g__Pseudomonas_E')


# In[13]:


plot_metaT_profiles('g__Methanothrix')


# In[15]:


def get_ratio_of_top_two(day):
    mean_day = jag.loc[(jag['treatment'] == 'catechin') & (jag['MAG'].str.startswith('contig')) & (jag['day'] == day)].groupby('MAG').agg({'getmm': 'mean'}).reset_index()
    
    # get ratio between average contig_591846 value and next highest vOTU
    ratio = mean_day.loc[mean_day['MAG'] == 'contig_591846']['getmm'].values[0] / mean_day.loc[mean_day['MAG'] != 'contig_591846']['getmm'].max()
    return ratio


# In[16]:


get_ratio_of_top_two(14), get_ratio_of_top_two(21), get_ratio_of_top_two(35)


# In[14]:


get_ipython().system('jupyter nbconvert --to script 007-vOTU_vs_host_metaT.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/04-A_vOTU_MAG_profiles_metaT')


# In[ ]:




