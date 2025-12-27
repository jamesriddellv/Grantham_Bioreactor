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


len(df)


# In[11]:


df['enzyme'].unique()


# In[12]:


# Get all MAGs that have either PGR or pgthAB (excluding those with FCR, PHY, CHI)
catechin_degraders = set(df.loc[df['enzyme'].isin(['FCR', 'PHY', 'CHI'])]['MAG'].unique())
pgr_mags = set(df.loc[df['enzyme'] == 'PGR']['MAG'].unique()) - catechin_degraders
pgthab_mags = set(df.loc[df['enzyme'] == 'pgthAB']['MAG'].unique()) - catechin_degraders

# Get all unique MAGs (union of both)
all_mags = pgr_mags | pgthab_mags


# In[13]:


MAG_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_GTDB_getmm.tsv', sep='\t')
MAG_metaT.head()


# In[14]:


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



# In[15]:


df['geTMM'].min()


# In[16]:


hydrogenase = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/hydrogenase_expression.csv')
hydrogenase


# In[17]:


enriched_mags = (sorted(enriched_3plus_timepoints) + 
               sorted(enriched_2_timepoints) + 
               sorted(enriched_1_timepoint)
                )
hydrogenase_enriched_only = hydrogenase.loc[hydrogenase['MAG'].isin(enriched_mags)]
len(enriched_mags), hydrogenase_enriched_only.MAG.nunique()


# In[ ]:


hydrogenase_enriched_only = hydrogenase_enriched_only.merge(df[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination']].drop_duplicates(), how='left')


# In[ ]:


hydrogenase_enriched_only = hydrogenase_enriched_only.rename(columns={'hydrogenase_group': 'description'})


# In[ ]:


h_melted = hydrogenase_enriched_only.melt(id_vars=['gene', 'description', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination'],
                              var_name='sample', value_name='geTMM')
h_melted


# In[ ]:


h_melted['enzyme'] = h_melted['description'].apply(lambda x: x.split('_')[0])


# In[ ]:


h_melted = h_melted.groupby(['gene', 'enzyme', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'sample', 'Completeness', 'Contamination']).agg({"geTMM": "sum"}).reset_index()


# In[ ]:


h_melted.columns


# In[ ]:


h_melted = h_melted.merge(df[['sample', 'treatment', 'time']].drop_duplicates(), on='sample', how='left')
h_melted.columns


# In[ ]:


df = df[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[ ]:


h_melted = h_melted[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[ ]:


df_merged = pd.concat([h_melted, df])


# In[ ]:


df_merged


# In[ ]:


# Convert to dictionary with MAG as key and (completeness, contamination) tuple as value
completeness_contamination_dict = dict(zip(
    completeness['MAG'],
    zip(completeness['Completeness'], completeness['Contamination'])
))


# In[ ]:


df_merged = df_merged.loc[(df_merged['Completeness'] > 50) & (df_merged['Contamination'] < 10)]
df_merged


# # Make heatmap with df_merged

# In[ ]:





# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Define enzymes to include
enzymes = ['PGR', 'pgthAB', 'carABCDE', 'CLD', 'mhpABCD', 'hpaAB', 'hpaDEFGH', 'fefe', 'nife']

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
    
    plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_enzyme_response_heatmap_enriched_only.pdf',
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


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 009-PGR-pgthAB-expression.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S9-PGR-pgthAB-expression')


# In[ ]:




