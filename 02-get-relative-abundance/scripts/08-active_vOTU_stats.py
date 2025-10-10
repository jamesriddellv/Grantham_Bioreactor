#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
vOTU_metadata.head()


# In[3]:


df = vOTU_metadata.drop_duplicates(subset='vOTU')


# In[4]:


# filter for only active vOTUs
getmms_vOTU_corrected_wide = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
active_vOTUs = list(getmms_vOTU_corrected_wide.loc[getmms_vOTU_corrected_wide.iloc[:,1:27].sum(axis=1) > 0].vOTU.unique())

df = df.loc[df['vOTU'].isin(active_vOTUs)]


# In[5]:


df.columns


# In[6]:


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


# In[7]:


df['completeness'] = df['completeness'].fillna(0)


# In[8]:


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


# In[9]:


df.describe()


# In[10]:


# check contig_591846
df.loc[df['vOTU'] == 'contig_591846']['completeness']


# In[11]:


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


# In[12]:


vOTU_metadata.loc[vOTU_metadata['vOTU'].isin(top_vOTUs)].sort_values(by=['completeness'])


# In[13]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_861080']


# In[14]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_776767']


# In[15]:


vOTU_metadata.loc[vOTU_metadata['vOTU']=='contig_773239']


# In[16]:


df=df.select_dtypes(exclude=['object'])


# In[17]:


df


# In[18]:


df.columns


# In[19]:


df = df[['length', 'n_genes', 'virus_score',
       'n_hallmarks', 'marker_enrichment',
       'viral_genes', 'host_genes', 'completeness',
       'contamination']]


# In[20]:


df


# In[21]:


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


# In[22]:


result = vOTU_metadata.groupby('vOTU').agg({
    'n_genes': 'first',  # or 'mean' since they're all the same per vOTU
    'concat_annotations': lambda x: x.str.contains('[a-zA-Z]', na=False).sum()
}).reset_index()

# Rename the column for clarity
result.columns = ['vOTU', 'n_genes', 'n_annotations_with_letters']

# Get the average number of genes
average_genes = result['n_genes'].mean()


# In[25]:


result['prop_annotated'] = result['n_annotations_with_letters'] / result['n_genes']
result['prop_annotated'].mean()


# In[26]:


result


# In[ ]:




