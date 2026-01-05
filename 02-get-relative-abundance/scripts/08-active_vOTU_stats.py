#!/usr/bin/env python
# coding: utf-8

# In[46]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[47]:


vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
vOTU_metadata.head()


# In[48]:


df = vOTU_metadata.drop_duplicates(subset='vOTU')


# In[49]:


# filter for only active vOTUs
getmms_vOTU_corrected_wide = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
active_vOTUs = list(getmms_vOTU_corrected_wide.loc[getmms_vOTU_corrected_wide.iloc[:,1:27].sum(axis=1) > 0].vOTU.unique())

df = df.loc[df['vOTU'].isin(active_vOTUs)]


# In[50]:


df.columns


# In[51]:


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Define vibrant, distinct colors
vibrant_colors = ['#FF6B35', '#004E89', '#00C9A7', '#845EC2', '#FFD23F']

# Create the boxplot with vibrant colors
plt.figure(figsize=(8,6))
ax = sns.boxplot(data=df, x='length', y='checkv_quality', 
                 # palette=vibrant_colors[:len(df['checkv_quality'].unique())])
                 color='white', linecolor='black', saturation=1)

# Get counts for each category and add to y-axis labels
category_counts = df['checkv_quality'].value_counts()
y_labels = []
for category in ax.get_yticklabels():
    category_name = category.get_text()
    count = category_counts.get(category_name, 0)
    y_labels.append(f"{category_name}\n(n={count})")

ax.set_yticklabels(y_labels)

# Remove legend
ax.legend().set_visible(False)

# Set x-axis to log10 scale
plt.xscale('log')
# Set specific ticks including zero
# ax.set_xticks([10000, 100000, 200000])

# Improve aesthetics
plt.xlabel('Active vOTU Length (log10 scale)', fontsize=16)
plt.ylabel('CheckV Quality', fontsize=16)
plt.xticks(fontsize=16)
# Increase tick length
ax.tick_params(axis='x', which='major', length=15)  # Default is usually 4
ax.tick_params(axis='x', which='minor', length=10)   # For minor ticks if you have them
plt.yticks(fontsize=16)

# Add mean line with label for legend
mean_line = plt.axvline(np.median(df['length']), color='red', label=f'Median Length 10.3kb')

# Add legend
plt.legend(handles=[mean_line], fontsize=10)

# plt.title('Distribution of Active vOTU Length by CheckV Quality', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1B_checkv_quality_vs_length.pdf', dpi=300)
plt.show()


# In[52]:


df['completeness'] = df['completeness'].fillna(0)


# In[53]:


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Use the same vibrant, distinct colors
vibrant_colors = ['#FF6B35', '#004E89', '#00C9A7', '#845EC2', 'grey']

# Get counts for each category
category_counts = df['checkv_quality'].value_counts()

# Create the scatterplot with vibrant colors
plt.figure(figsize=(12, 8))
ax = sns.scatterplot(data=df, x='length', y='completeness', 
                     hue='checkv_quality',
                     palette=vibrant_colors[:len(df['checkv_quality'].unique())],
                     alpha=0.7, s=60)

# Set x-axis to log10 scale
plt.xscale('log')

# Increase tick label sizes
ax.tick_params(axis='both', which='major', labelsize=16)

# Get the legend and modify labels to include counts
legend = ax.legend(fontsize=16, title_fontsize=18, markerscale=2)
legend.set_title('CheckV Quality', prop={'size': 18, 'weight': 'bold'})

# Update legend labels to include counts
for text in legend.get_texts():
    category = text.get_text()
    count = category_counts[category]
    text.set_text(f'{category} (n={count})')

# Improve aesthetics
plt.xlabel('vOTU Length', fontsize=18)
plt.ylabel('Checkv Completeness', fontsize=18)
plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-C-completeness-vs-length.pdf', dpi=300)
plt.show()


# In[54]:


df.describe()


# In[55]:


# check contig_591846
df.loc[df['vOTU'] == 'contig_591846']['completeness']


# In[56]:


top_vOTUs = ['contig_1022629',
 'contig_1093732',
 'contig_507316',
 'contig_588865',
 'contig_591846',
 'contig_709158',
 'contig_713418',
 'contig_724521',
 'contig_755337',
 'contig_773239',
 'contig_776767',
 'contig_861080']


# In[57]:


top_vOTUs = [
    'contig_1022629',
    'contig_1093732',
    'contig_709158',
    'contig_713418',
    'contig_773239',
    'contig_776767',
    'contig_861080'
]


# In[58]:


colors_dict_fig_2 = {
'no_assignment__contig_1022629': '#C0C0C0',     # Silver (from original)
    'no_assignment__contig_1093732': '#D3D3D3',     # Light gray (from original)
    'no_assignment__contig_709158': '#696969',      # Dim gray (from original)
    'no_assignment__contig_713418': '#778899',      # Light slate gray (from original)
    'no_assignment__contig_773239': '#708090',      # Slate gray (from original)
    'no_assignment__contig_776767': '#556B2F',      # Dark olive green (from original)
    'no_assignment__contig_861080': '#2F4F4F',      # Dark slate gray (from original)
}


# In[59]:


vOTU_metadata.loc[vOTU_metadata['vOTU'].isin(top_vOTUs)].sort_values(by=['completeness'])


# In[60]:


vOTU_metadata.loc[vOTU_metadata['vOTU'].isin(top_vOTUs)].sort_values(by=['completeness']).drop_duplicates(subset='vOTU')


# In[61]:


vOTU_metadata.loc[vOTU_metadata['vOTU'].isin(top_vOTUs)].sort_values(by=['completeness']).drop_duplicates(subset='vOTU').describe()


# In[62]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_861080']


# In[63]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_776767']


# In[64]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_773239']


# In[65]:


df=df.select_dtypes(exclude=['object'])


# In[66]:


df


# In[67]:


df.columns


# In[68]:


df = df[['length', 'n_genes', 'virus_score',
       'n_hallmarks', 'marker_enrichment',
       'viral_genes', 'host_genes', 'completeness',
       'contamination']]


# In[69]:


df


# In[70]:


# Create a much larger figure
plt.figure(figsize=(30, 30))

# Create scatter matrix with adjusted parameters
scatter_matrix = pd.plotting.scatter_matrix(
    df, 
    alpha=0.2,
    figsize=(30, 30),
    diagonal='hist',  # or 'kde' for kernel density estimation
    hist_kwds={'bins': 50}
)

# Adjust font sizes for better readability
for ax in scatter_matrix.flatten():
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontsize(12)
    ax.tick_params(labelsize=10)

# Adjust layout to prevent overlapping
plt.tight_layout(pad=2.0)

# Save with high DPI
# plt.savefig('scatter_matrix_active_vOTUs_large.pdf', dpi=300, bbox_inches='tight')
plt.show()


# In[71]:


result = vOTU_metadata.groupby('vOTU').agg({
    'n_genes': 'first',  # or 'mean' since they're all the same per vOTU
    'concat_annotations': lambda x: x.str.contains('[a-zA-Z]', na=False).sum()
}).reset_index()

# Rename the column for clarity
result.columns = ['vOTU', 'n_genes', 'n_annotations_with_letters']

# Get the average number of genes
average_genes = result['n_genes'].mean()


# In[72]:


result['prop_annotated'] = result['n_annotations_with_letters'] / result['n_genes']
result['prop_annotated'].mean()


# # For the unassigned vOTUs, are they related to any vOTUs with a host prediction?

# In[73]:


no_assignment_vOTUs = [
    'contig_1022629',
    'contig_1093732',
    'contig_709158',
    'contig_713418',
    'contig_773239',
    'contig_776767',
    'contig_861080'
]


# In[74]:


# add iphop and other metadata
iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df.head()


# In[77]:


vcontact3_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-3.1.7-out/exports/final_assignments_with_pipe_chars.csv')
vcontact3_df = vcontact3_df.fillna('')
vcontact3_df.head()


# In[78]:


related_vOTUs = {}

for i in no_assignment_vOTUs:

    
    
    genus = list(vcontact3_df.loc[vcontact3_df['genus_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['genus_prediction'].values[0]]['Genome'])
    subfamily = list(vcontact3_df.loc[vcontact3_df['subfamily_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['subfamily_prediction'].values[0]]['Genome'])
    family = list(vcontact3_df.loc[vcontact3_df['family_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['family_prediction'].values[0]]['Genome'])

    if vcontact3_df.loc[vcontact3_df['Genome'] == i]['realm_prediction'].values[0] == 'singleton':
        related_vOTUs[i] = ['singleton']
    
    elif len(genus) > 1:
        related_vOTUs[i] = ['genus'] + genus
        
    elif len(subfamily) > 1:
        related_vOTUs[i] = ['subfamily'] + subfamily
        
    elif len(family) > 1:
        related_vOTUs[i] = ['family'] + family
    else:
        related_vOTUs[i] = ['no_relatives']

related_vOTUs


# In[79]:


iphop_df.loc[iphop_df['vOTU'].isin(['contig_1004352|provirus_3_15586',
  'contig_675166',
  'contig_764965',
  'contig_773239'])]


# In[80]:


# now for each of these values, get the iphop results:

related_vOTU_predictions = {}

for k,v in related_vOTUs.items():
    related_vOTU_predictions[k] = set(iphop_df.loc[iphop_df['vOTU'].isin(v)]['highest_host_tax_rank'])

related_vOTU_predictions


# In[81]:


iphop_df.loc[iphop_df['vOTU'].isin(no_assignment_vOTUs)]


# In[82]:


related_vOTUs = {}

for i in no_assignment_vOTUs:
    genus = list(vcontact3_df.loc[vcontact3_df['genus_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['genus_prediction'].values[0]]['Genome'])
    subfamily = list(vcontact3_df.loc[vcontact3_df['subfamily_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['subfamily_prediction'].values[0]]['Genome'])
    family = list(vcontact3_df.loc[vcontact3_df['family_prediction'] == vcontact3_df.loc[vcontact3_df['Genome'] == i]['family_prediction'].values[0]]['Genome'])

    if vcontact3_df.loc[vcontact3_df['Genome'] == i]['realm_prediction'].values[0] == 'singleton':
        related_vOTUs[i] = ['singleton']
    
    elif len(genus) > 1:
        related_vOTUs[i] = genus
        
    elif len(subfamily) > 1:
        related_vOTUs[i] = subfamily
        
    elif len(family) > 1:
        related_vOTUs[i] = family

related_vOTUs


# In[83]:


unassigned_vOTU_relatives = [item for sublist in related_vOTUs.values() for item in sublist]


# In[84]:


unassigned_vOTU_relatives_tax = vcontact3_df.loc[vcontact3_df['Genome'].isin(unassigned_vOTU_relatives)]


# In[85]:


unassigned_vOTU_relatives_tax = unassigned_vOTU_relatives_tax.rename(columns={'Genome': 'vOTU'})


# In[86]:


unassigned_vOTU_relatives_tax.columns


# In[87]:


vcontact3_df.loc[vcontact3_df['Genome'] == 'contig_861080']


# In[88]:


unassigned_vOTU_relatives_taxa_predicted_hosts = (
    unassigned_vOTU_relatives_tax
    .merge(iphop_df, on='vOTU', how='left')
    .fillna('no_assignment')
    .drop_duplicates(subset=['vOTU', 'highest_host_tax_rank'])
    .groupby(['vOTU', 
       'realm_prediction', 
       'kingdom_prediction', 'phylum_prediction',
       'class_prediction',
       'order_prediction', 'family_prediction',
       'subfamily_prediction',
       'genus_prediction', 'MAG_is_active', 'MAG_is_stordalen'])
    .agg({'highest_host_tax_rank': ';'.join})
).reset_index()

unassigned_vOTU_relatives_taxa_predicted_hosts = unassigned_vOTU_relatives_taxa_predicted_hosts.merge(vOTU_metadata[['vOTU', 'length', 'n_genes', 'completeness']].drop_duplicates(), on='vOTU', how='left')
unassigned_vOTU_relatives_taxa_predicted_hosts['length (kb)'] = unassigned_vOTU_relatives_taxa_predicted_hosts['length'] / 1000
unassigned_vOTU_relatives_taxa_predicted_hosts


# In[89]:


unassigned_vOTU_relatives_taxa_predicted_hosts['vOTU'] = unassigned_vOTU_relatives_taxa_predicted_hosts['vOTU'].str.replace('|', '_')
unassigned_vOTU_relatives_taxa_predicted_hosts['vOTU']


# In[90]:


unassigned_vOTU_relatives_taxa_predicted_hosts.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/unassigned_vOTU_relatives.tsv', sep='\t', index=False)


# In[91]:


# load vcontact3 network and subset these files
edge_table = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-out/exports/vcontact3_edge_table.csv')
edge_table = edge_table.loc[
(edge_table['source'].isin(unassigned_vOTU_relatives_taxa_predicted_hosts['vOTU']))
 & (edge_table['target'].isin(unassigned_vOTU_relatives_taxa_predicted_hosts['vOTU']))
 ]
edge_table.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/unassigned_vOTU_relatives_vcontact3_edge_table.csv', index=False)


# In[25]:


get_ipython().system('jupyter nbconvert --to script 006-active_vOTU_stats.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts/08-active_vOTU_stats')


# In[ ]:




