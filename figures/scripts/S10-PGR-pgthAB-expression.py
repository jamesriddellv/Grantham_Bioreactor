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


# # Now figure out if these pathways were active, and if they were upregulated in catechin-amended samples.

# In[17]:


# load per-gene activity for vOTUs, across all samples, isolate these genes.
genes_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/geTMMs_ge5_1X_REV_26samples.txt', sep='\t')
genes_df


# In[18]:


camper_genes_df = genes_df.loc[genes_df['gene'].isin(active_camper_annotations['gene'].unique())]
camper_genes_df


# In[19]:


camper_genes_df = camper_genes_df.merge(active_camper_annotations, on='gene', how='left')


# In[20]:


camper_genes_df['MAG'] = camper_genes_df['fasta']


# In[21]:


# Get all MAGs that have either PGR or pgthAB (excluding those with FCR, PHY, CHI)
catechin_degraders = set(df.loc[df['enzyme'].isin(['FCR', 'PHY', 'CHI'])]['MAG'].unique())
pgr_mags = set(df.loc[df['enzyme'] == 'PGR']['MAG'].unique()) - catechin_degraders
pgthab_mags = set(df.loc[df['enzyme'] == 'pgthAB']['MAG'].unique()) - catechin_degraders

# Get all unique MAGs (union of both)
all_mags = pgr_mags | pgthab_mags


# In[22]:


MAG_metaT = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_GTDB_getmm.tsv', sep='\t')
MAG_metaT.head()


# In[23]:


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



# In[24]:


df['geTMM'].min()


# In[25]:


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

# In[26]:


hydrogenase = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/hydrogenase_expression.csv')
hydrogenase


# In[27]:


enriched_mags = (sorted(enriched_3plus_timepoints) + 
               sorted(enriched_2_timepoints) + 
               sorted(enriched_1_timepoint)
                )
hydrogenase_enriched_only = hydrogenase.loc[hydrogenase['MAG'].isin(enriched_mags)]
len(enriched_mags), hydrogenase_enriched_only.MAG.nunique()


# In[28]:


hydrogenase_enriched_only = hydrogenase_enriched_only.merge(df[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination']].drop_duplicates(), how='left')


# In[29]:


hydrogenase_enriched_only = hydrogenase_enriched_only.rename(columns={'hydrogenase_group': 'description'})


# In[30]:


h_melted = hydrogenase_enriched_only.melt(id_vars=['gene', 'description', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination'],
                              var_name='sample', value_name='geTMM')
h_melted


# In[31]:


h_melted['enzyme'] = h_melted['description'].apply(lambda x: x.split('_')[0])


# In[32]:


h_melted = h_melted.groupby(['gene', 'enzyme', 'MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'sample', 'Completeness', 'Contamination']).agg({"geTMM": "sum"}).reset_index()


# In[33]:


h_melted.columns


# In[34]:


h_melted = h_melted.merge(df[['sample', 'treatment', 'time']].drop_duplicates(), on='sample', how='left')
h_melted.columns


# In[35]:


df = df[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[36]:


h_melted = h_melted[['gene', 'MAG', 'enzyme', 'sample',
       'treatment', 'time', 'geTMM', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S',
       'highest_host_tax_rank', 'Completeness', 'Contamination']]


# In[37]:


h_melted.MAG.nunique()


# In[38]:


df_merged = pd.concat([h_melted, df])


# In[39]:


df_merged


# In[40]:


# Convert to dictionary with MAG as key and (completeness, contamination) tuple as value
completeness_contamination_dict = dict(zip(
    completeness['MAG'],
    zip(completeness['Completeness'], completeness['Contamination'])
))


# In[41]:


df_merged = df_merged.loc[(df_merged['Completeness'] > 50) & (df_merged['Contamination'] < 10)]
df_merged


# # Add additional polpyhenol metabolisms

# In[42]:


df_merged.columns


# In[43]:


camper_genes_df = camper_genes_df.merge(df[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank', 'Completeness', 'Contamination']].drop_duplicates(), how='left', on='MAG')
camper_genes_df = camper_genes_df.loc[camper_genes_df['MAG'].isin(enriched_mags)]
camper_genes_df


# In[44]:


camper_genes_df.columns


# In[45]:


camper_genes_melted = camper_genes_df.melt(id_vars=['gene', 'fasta', 'scaffold',
       'gene_position', 'camper_hits', 'camper_rank', 'camper_bitScore',
       'camper_id', 'camper_definition', 'camper_search_type', 'MAG', 'GTDB',
       'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank',
       'Completeness', 'Contamination'], var_name='sample', value_name='geTMM')

camper_genes_melted['enzyme'] = camper_genes_melted['camper_definition'].fillna('').apply(lambda x: x.split(';', 1)[0])


# In[46]:


camper_genes_melted


# In[47]:


camper_genes_melted.enzyme.unique()


# In[48]:


camper_genes_melted.camper_definition.unique()


# In[49]:


# filter for pgthAB and PGR-containing MAGs
# pgthab_pgr_mags = pgr_mags | pgthab_mags
# pgthab_pgr_mags_genes_melted = camper_genes_melted.loc[camper_genes_melted['fasta'].isin(pgthab_pgr_mags)]


# In[50]:


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


# In[51]:


additional_melted = camper_genes_melted.loc[camper_genes_melted['enzyme'].isin(additional_enzymes)]
additional_melted


# In[52]:


additional_melted = additional_melted.merge(df_merged[['sample', 'treatment', 'time']].drop_duplicates(), on='sample', how='left')


# In[53]:


additional_melted.MAG.nunique()


# In[54]:


df_merged.columns


# In[55]:


additional_melted_to_merge = additional_melted[['gene', 'MAG', 'enzyme', 'sample', 'treatment', 'time', 'geTMM', 'GTDB',
       'K', 'P', 'C', 'O', 'F', 'G', 'S', 'highest_host_tax_rank',
       'Completeness', 'Contamination']]


# In[56]:


df_merged.MAG.nunique()


# In[57]:


df_merged = pd.concat([df_merged, additional_melted_to_merge])


# # Make heatmap with df_merged

# In[58]:


import os
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

def plot_foldchange_heatmap_enriched_only(fc_df, mag_to_tax, export_dir=None):
    """
    Create heatmap of log2 fold-changes with timepoints for each enzyme.
    Enzymes on y-axis, MAGs on x-axis, with timepoints as sub-columns.
    Additionally, converts matrix values to a long-form tidy table for compliance.
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
    
    # -----------------------------------------------------------------
    # DATA POLICY EXPORT BLOCK (Generates tidy, peer-review-ready CSV)
    # -----------------------------------------------------------------
    if export_dir:
        os.makedirs(export_dir, exist_ok=True)
        source_data_records = []
        
        # Unroll the matrix into tidy coordinate pairs matching the figure grid
        for mag in column_order:
            mag_tax = mag_to_tax.get(mag, 'Unknown')
            mag_stats = completeness_contamination_dict.get(mag, 'Unknown')
            full_resolved_label = f"{mag}; {mag_tax}; {mag_stats}"
            group_type = fc_df.loc[mag, 'Enrichment_Type']
            
            for row_idx in heatmap_data.index:
                enzyme_name = row_idx.split('_')[0]
                day_val = row_idx.split('_day')[1]
                log2_fc_val = heatmap_data.loc[row_idx, mag]
                
                source_data_records.append({
                    'MAG_ID': mag,
                    'Taxonomy_and_Stats': full_resolved_label,
                    'Enrichment_Cohort': f"Catechin-enriched ({group_type})",
                    'Enzyme': enzyme_name,
                    'Day': int(day_val),
                    'Log2_Fold_Change': log2_fc_val if not pd.isna(log2_fc_val) else "NA"
                })
        
        df_export = pd.DataFrame(source_data_records)
        export_file_path = os.path.join(export_dir, '5A_MAG_polyphenol_activity_heatmap_source_data.csv')
        df_export.to_csv(export_file_path, index=False)
        print(f"\n[Data Policy] Tidy Source Data written successfully to: {export_file_path}")
    # -----------------------------------------------------------------

    # Create figure
    fig = plt.figure(figsize=(max(16, len(column_order) * 0.5), 10))
    ax = fig.add_subplot(111)
    
    from matplotlib.colors import LinearSegmentedColormap
    colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
              '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
    cmap.set_bad(color='#e0e0e0')  # Light gray for missing data
    
    # Determine vmin and vmax for symmetric color scale
    vmax = np.nanmax(np.abs(heatmap_data.values))
    vmin = -vmax
    
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
    
    for i in range(len(enzymes)):
        if i > 0:
            ax.axhline(y=i*len(timepoints), color='black', linewidth=3)
    
    for i in range(len(column_order)):
        if i > 0:
            ax.axvline(x=i, color='black', linewidth=1.5)
    
    ax.set_ylabel('Enzyme_Day', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9, length=0)
    
    ax.tick_params(axis='x', labelsize=7, rotation=45, labelbottom=False, labeltop=True, 
                    bottom=False, top=True)
    ax.xaxis.set_label_position('top')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    
    boundaries = []
    prev_type = None
    for i, enrich_type in enumerate(enrichment_types):
        if enrich_type != prev_type and prev_type is not None:
            boundaries.append(i)
        prev_type = enrich_type
    
    for boundary in boundaries:
        ax.axvline(x=boundary, color='white', linewidth=4, linestyle='-')
    
    group_starts = [0] + boundaries
    group_ends = boundaries + [len(column_order)]
    
    for start, end in zip(group_starts, group_ends):
        midpoint = (start + end) / 2
        enrich_type = enrichment_types[start]
        ax.text(midpoint, len(heatmap_data) + 1.2, f'Catechin-enriched\n{enrich_type}',
               va='bottom', ha='center', fontsize=10, fontweight='bold',
               rotation=0)
    
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

# --- Run pipeline with direct export variable ---
target_data_directory = '/users/PAS1573/riddell26/data/'

fc_df_enriched = calculate_log2_foldchange_heatmap_enriched_only(
    enriched_3plus_timepoints, 
    enriched_2_timepoints, 
    enriched_1_timepoint,
    mag_to_tax, 
    df_merged
)

fig = plot_foldchange_heatmap_enriched_only(fc_df_enriched, mag_to_tax, export_dir=target_data_directory)


# In[59]:


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


# In[60]:


# Define enzymes to include
enzymes = ['FLR', 'CHI', 'FCR', 'PHY', 'PGR', 'pgthAB', # Flavonoids
           'fefe', 'nife', # Hydrogenases
           'CLD', 'EDL', 'PhcD', 'BER', 'SesA', 'PinZ', 'PAH', 'GLM', 'PhcC', 'PhcG', # Lignans
           'EhyB', 'CalA', 'CalB', 'IEM', # Phenyl-propenes
           'FdeE', 'FdeC', # Flavonones
           'CurA', 'carABCDE', 'mhpABCD', 'hpaAB', 'hpaDEFGH', 'GRAD', # Phenolic acids
           'sam5', # Isoflavonoids, cinnamic acids
           'ChlE', # Cinnamic acid metabolism
           'FCS'] # Non-specific

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

def plot_foldchange_heatmap_enriched_only(fc_df, mag_to_tax, export_dir=None):
    """
    Create heatmap of log2 fold-changes with timepoints for each enzyme.
    Enzymes on y-axis, MAGs on x-axis, with timepoints as sub-columns.
    Additionally, converts matrix values to a long-form tidy table for compliance.
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
    
    # -----------------------------------------------------------------
    # DATA POLICY EXPORT BLOCK (Generates S13_ tidy source data table)
    # -----------------------------------------------------------------
    if export_dir:
        os.makedirs(export_dir, exist_ok=True)
        source_data_records = []
        
        # Unroll the matrix into tidy coordinate pairs matching the expanded dataset
        for mag in column_order:
            mag_tax = mag_to_tax.get(mag, 'Unknown')
            mag_stats = completeness_contamination_dict.get(mag, 'Unknown')
            full_resolved_label = f"{mag}; {mag_tax}; {mag_stats}"
            group_type = fc_df.loc[mag, 'Enrichment_Type']
            
            for row_idx in heatmap_data.index:
                enzyme_name = row_idx.split('_')[0]
                day_val = row_idx.split('_day')[1]
                log2_fc_val = heatmap_data.loc[row_idx, mag]
                
                source_data_records.append({
                    'MAG_ID': mag,
                    'Taxonomy_and_Stats': full_resolved_label,
                    'Enrichment_Cohort': f"Catechin-enriched ({group_type})",
                    'Enzyme': enzyme_name,
                    'Day': int(day_val),
                    'Log2_Fold_Change': log2_fc_val if not pd.isna(log2_fc_val) else "NA"
                })
        
        df_export = pd.DataFrame(source_data_records)
        export_file_path = os.path.join(export_dir, 'S13_additional_metabolisms_heatmap_source_data.csv')
        df_export.to_csv(export_file_path, index=False)
        print(f"\n[Data Policy] Tidy Source Data written successfully to: {export_file_path}")
    # -----------------------------------------------------------------

    # Create figure
    fig = plt.figure(figsize=(max(16, len(column_order) * 0.5), 20))
    ax = fig.add_subplot(111)
    
    from matplotlib.colors import LinearSegmentedColormap
    colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
              '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
    cmap.set_bad(color='#e0e0e0')  # Light gray for missing data
    
    # Determine vmin and vmax for symmetric color scale
    vmax = np.nanmax(np.abs(heatmap_data.values))
    vmin = -vmax
    
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
    
    for i in range(len(enzymes)):
        if i > 0:
            ax.axhline(y=i*len(timepoints), color='black', linewidth=3)
    
    for i in range(len(column_order)):
        if i > 0:
            ax.axvline(x=i, color='black', linewidth=1.5)
    
    ax.set_ylabel('Enzyme_Day', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=9, length=0)
    
    ax.tick_params(axis='x', labelsize=7, rotation=45, labelbottom=False, labeltop=True, 
                    bottom=False, top=True)
    ax.xaxis.set_label_position('top')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='left', rotation_mode='anchor')
    
    boundaries = []
    prev_type = None
    for i, enrich_type in enumerate(enrichment_types):
        if enrich_type != prev_type and prev_type is not None:
            boundaries.append(i)
        prev_type = enrich_type
    
    for boundary in boundaries:
        ax.axvline(x=boundary, color='white', linewidth=4, linestyle='-')
    
    group_starts = [0] + boundaries
    group_ends = boundaries + [len(column_order)]
    
    for start, end in zip(group_starts, group_ends):
        midpoint = (start + end) / 2
        enrich_type = enrichment_types[start]
        ax.text(midpoint, len(heatmap_data) + 1.2, f'Catechin-enriched\n{enrich_type}',
               va='bottom', ha='center', fontsize=10, fontweight='bold',
               rotation=0)
    
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

# --- Run pipeline with direct export variable ---
target_data_directory = '/users/PAS1573/riddell26/data/'

fc_df_enriched = calculate_log2_foldchange_heatmap_enriched_only(
    enriched_3plus_timepoints, 
    enriched_2_timepoints, 
    enriched_1_timepoint,
    mag_to_tax, 
    df_merged
)

fig = plot_foldchange_heatmap_enriched_only(fc_df_enriched, mag_to_tax, export_dir=target_data_directory)


# In[61]:


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


# # Investigate Gallate -> pyrogallol metabolism

# In[62]:


getmm_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/geTMM_table.csv')
getmm_df.head()


# In[63]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[64]:


lpdC_annotations = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/gallate_decarboxylase_annotations.tsv', sep='\t')
lpdC_annotations = lpdC_annotations.rename(columns={'Unnamed: 0': 'gene'})
lpdC_annotations


# In[65]:


lpdC_annotations['pfam_hits'].value_counts()


# In[66]:


# Of the 52 genes encoding lpdC, 20 were also identified in the getmm table.
lpdC_activity = getmm_df.loc[getmm_df.gene.isin(lpdC_annotations.gene)]
lpdC_activity


# In[67]:


non_cat_pgth_degraders = lpdC_activity.loc[lpdC_activity.fasta.isin(pgthab_mags)]

non_cat_pgth_degraders.drop(columns='fasta', inplace=True)


df_melted = non_cat_pgth_degraders.melt(id_vars='gene', var_name='Sample', value_name='getmm')
df_melted = df_melted.merge(replicate_frame, on='Sample', how='left')
df_melted


# In[68]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Select four genes to plot
genes_to_plot = df_melted['gene'].unique()  # Adjust this to select specific genes

# Create a 2x2 subplot layout
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

# Define colors for treatments
colors = {'catechin': '#FC9B2D', 'unamended': '#7ACBC3'}

# Plot each gene
for idx, gene in enumerate(genes_to_plot):
    ax = axes[idx]
    
    # Filter data for this gene
    gene_data = df_melted[df_melted['gene'] == gene]
    
    # Calculate mean and std for each treatment-day combination
    summary = gene_data.groupby(['treatment', 'day'])['getmm'].agg([
        ('mean', 'mean'),
        ('std', 'std'),
        ('count', 'count')
    ]).reset_index()
    
    # Plot each treatment
    for treatment in summary['treatment'].unique():
        treatment_data = summary[summary['treatment'] == treatment]
        
        # Plot replicate points
        replicate_data = gene_data[gene_data['treatment'] == treatment]
        ax.scatter(replicate_data['day'], replicate_data['getmm'],
                  color=colors.get(treatment, 'gray'), alpha=0.6, s=50, zorder=3)
        
        # Plot mean line
        ax.plot(treatment_data['day'], treatment_data['mean'], 
                linewidth=2.5,
                label=treatment, color=colors.get(treatment, 'gray'))
    
    # Formatting
    ax.set_xlabel('Day', fontsize=11, fontweight='bold')
    ax.set_ylabel('getmm', fontsize=11, fontweight='bold')
    ax.set_title(f'{gene}', 
                 fontsize=9, fontweight='bold')
    ax.legend(title='Treatment', fontsize=9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks([0, 7, 14, 21, 35])

plt.tight_layout()
plt.show()


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 009-S10-PGR-pgthAB-expression.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S10-PGR-pgthAB-expression')

