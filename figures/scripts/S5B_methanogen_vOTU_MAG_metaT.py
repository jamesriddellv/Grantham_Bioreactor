#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
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


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df


# In[6]:


methanogen_vOTUs = iphop_df.loc[iphop_df['highest_host_tax_rank'].str.contains('g__Methano.*', na=False)]['vOTU'].unique()


# In[7]:


vOTU_metaT = vOTU_metaT.drop(columns=['highest_host_tax_rank']).merge(iphop_df[['vOTU', 'highest_host_tax_rank']].drop_duplicates(subset='vOTU'), on='vOTU', how='left')
vOTU_metaT


# In[8]:


methanogen_vOTUs = vOTU_metaT.loc[vOTU_metaT['highest_host_tax_rank'].str.contains('g__Methano.*', na=False)]
methanogen_vOTUs


# In[9]:


color_dict = {
    'contig_1091961_g__Methanosarcina': '#8B5FBF',  # Deep purple
    'contig_1166240_g__Methanothrix': '#00CED1',    # Dark turquoise
    'contig_391099_g__Methanobacterium_A': '#FF6B6B', # Coral red
    'contig_535315_g__Methanothrix': '#4ECDC4',     # Mint teal
    'contig_771149_g__Methanobacterium_A': '#FFA62B', # Bright orange
    'contig_933694_g__Methanosarcina': '#6A0572',   # Deep magenta
}



# In[12]:


def plot_metaT_profiles_by_taxa(host):

    # filter to host only
    host_df = MAG_metaT.loc[MAG_metaT['highest_host_tax_rank'].str.contains('g__Methan.*', na=False)]
    host_df_melted = host_df.iloc[:,:27].melt(id_vars='MAG', var_name='Sample', value_name='getmm')
    # Merge back the taxonomy info
    host_df_melted = host_df_melted.merge(
        host_df[['MAG', 'highest_host_tax_rank']], 
        on='MAG', 
        how='left'
    )

    vOTU_df = vOTU_metaT.loc[vOTU_metaT['highest_host_tax_rank'].str.contains('g__Methan.*', na=False)]
    vOTU_df_melted = vOTU_df.iloc[:,:27].melt(id_vars='vOTU', var_name='Sample', value_name='getmm')
    # Merge back the taxonomy info
    vOTU_df_melted = vOTU_df_melted.merge(
        vOTU_df[['vOTU', 'highest_host_tax_rank']], 
        on='vOTU', 
        how='left'
    )
    
    # load metadata
    metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv',
                      usecols=['Sample', 'treatment', 'timepoint'])
    metadata['day'] = metadata['timepoint'].apply(lambda x: int(x.replace('day', '')))
    metadata = metadata.drop(columns=['timepoint'])

    # merge to host and vOTU dfs
    host_df_melted = host_df_melted.merge(metadata, on='Sample', how='left')
    vOTU_df_melted = vOTU_df_melted.merge(metadata, on='Sample', how='left')

    # rename to match
    vOTU_df_melted.rename(columns={'vOTU': 'MAG'}, inplace=True)
    
    # order shared columns the same
    host_df_melted = host_df_melted[['MAG', 'Sample', 'treatment', 'day', 'getmm', 'highest_host_tax_rank']]
    vOTU_df_melted = vOTU_df_melted[['MAG', 'Sample', 'treatment', 'day', 'getmm', 'highest_host_tax_rank']]

    # merge dataframes
    merged_df = pd.concat([host_df_melted, vOTU_df_melted], axis=0)

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
    merged_df_mean = merged_df.groupby(['MAG', 'treatment', 'day', 'highest_host_tax_rank']).agg({'getmm': 'mean'}).reset_index()

    # add provirus and MAG designations for labeling
    merged_df_mean['is_provirus'] = merged_df_mean['MAG'].str.contains('provirus')
    merged_df_mean['is_MAG'] = ~merged_df_mean['MAG'].str.contains('contig_')

    # log transform GeTMM since values are a broad range
    merged_df_mean['log10p_getmm'] = merged_df_mean['getmm'].apply(lambda x: np.log10(x+1))

    merged_df_mean.to_csv('data/S5B_mean_methanogen_vOTU_MAG_getmm.tsv')
    
    # Get unique taxonomic ranks
    unique_taxa = sorted(merged_df_mean['highest_host_tax_rank'].unique())
    n_taxa = len(unique_taxa)
    
    # Assign a base color to each taxonomy
    base_colors = plt.cm.tab10(np.linspace(0, 1, n_taxa))
    taxa_base_colors = dict(zip(unique_taxa, base_colors))

    # Set style and figure size
    sns.set_style('whitegrid')
    fig, axes = plt.subplots(
        n_taxa, 2, figsize=(14, 4*n_taxa), sharex=True
    )  # n_taxa rows, 2 columns (unamended, catechin)
    
    # Ensure axes is 2D even if only one taxonomy
    if n_taxa == 1:
        axes = axes.reshape(1, -1)
    
    treatments = ['unamended', 'catechin']
    
    for row_idx, taxonomy in enumerate(unique_taxa):
        # Filter data for this taxonomy
        taxa_data = merged_df_mean[merged_df_mean['highest_host_tax_rank'] == taxonomy]
        
        # Get all unique MAGs for this taxonomy (sorted in reverse for legend)
        mags_in_taxa = sorted(taxa_data['MAG'].unique(), reverse=True)
        n_mags = len(mags_in_taxa)
        
        # Create color shades for each MAG using the base color
        base_color = taxa_base_colors[taxonomy]
        # Generate shades from lighter to darker
        from matplotlib.colors import to_rgb
        rgb = to_rgb(base_color)
        
        # Create shades by adjusting brightness (reverse order for color assignment)
        shades = []
        for i in range(n_mags):
            # Interpolate between lighter (0.3 offset) and darker (1.0) versions
            factor = 0.3 + (0.7 * i / max(n_mags - 1, 1))
            shade = tuple(c * factor for c in rgb)
            shades.append(shade)
        
        # Reverse mags_in_taxa for color assignment so first (bottom in legend) gets darkest
        mag_color_map = dict(zip(reversed(mags_in_taxa), shades))
        
        # Store handles and labels for legend
        legend_handles = []
        legend_labels = []
        
        for col_idx, treatment in enumerate(treatments):
            ax = axes[row_idx, col_idx]
            
            # Make a twin axis for MAGs
            ax_mag = ax.twinx()
            
            # Filter data for this treatment
            treatment_data = taxa_data[taxa_data['treatment'] == treatment]
            
            # Plot in reverse order so legend order matches visual order
            for mag_name in reversed(mags_in_taxa):
                # Filter data for this specific MAG and sort by day
                mag_data = treatment_data[treatment_data['MAG'] == mag_name].sort_values('day')
                
                if len(mag_data) == 0:
                    continue
                
                # First row for style decisions
                first_row = mag_data.iloc[0]
                
                linestyle = '--' if first_row['is_MAG'] else '-'
                linewidth = 3 if first_row['is_provirus'] or first_row['is_MAG'] else 1
                color = mag_color_map[mag_name]
                
                # Choose which axis to plot on
                target_ax = ax_mag if first_row['is_MAG'] else ax
                
                line, = target_ax.plot(
                    mag_data['day'],
                    mag_data['log10p_getmm'],
                    label=mag_name,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    marker='o',
                    markersize=4,
                    alpha=0.8,
                    color=color
                )
                
                # Collect handles for legend (only from catechin plot)
                if col_idx == 1 and mag_name not in legend_labels:
                    legend_handles.append(line)
                    legend_labels.append(mag_name)
            
            # Format axes
            ax.set_title(f"{taxonomy} - {treatment}", fontsize=11, fontweight='bold')
            ax.set_xticks([0, 14, 21, 35])
            ax.tick_params(axis='both', which='major', labelsize=10)
            
            # Only label y-axis on leftmost column
            if col_idx == 0:
                ax.set_ylabel('GeTMM (log10(x+1))', fontsize=11)
                ax_mag.tick_params(axis='y', colors="gray")
            
            # Label MAG axis on rightmost column
            if col_idx == 1:
                ax_mag.set_ylabel("MAGs (GeTMM log(x+1))", fontsize=11)
                ax_mag.tick_params(axis='y', colors="gray")
                
                # Add legend to the right of catechin plot
                ax_mag.legend(legend_handles, legend_labels, 
                             loc='center left', 
                             bbox_to_anchor=(1.15, 0.5),
                             fontsize=9,
                             frameon=True)
            else:
                ax_mag.tick_params(axis='y', colors="gray")
    
    # Shared x-axis label on bottom row
    for col_idx in range(2):
        axes[-1, col_idx].set_xlabel("Day", fontsize=12)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)  # Make room for legends
    plt.savefig(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4_methanogens_vOTU_metaT.pdf', dpi=300, bbox_inches='tight')
    plt.show()




# In[13]:


plot_metaT_profiles_by_taxa('methanogens')


# In[ ]:




