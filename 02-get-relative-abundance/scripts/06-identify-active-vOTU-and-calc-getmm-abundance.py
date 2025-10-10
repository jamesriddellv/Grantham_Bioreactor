#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# ## 1. Extract gene IDs from vOTU database that have a structural phage or lysis gene annotation

# In[2]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
df


# In[3]:


# Filter for virus structural or lysis genes
lytic_keywords = 'capsid|tail|packag|portal|lysis|terminase|spike|baseplate|internal virion|neck|prohead'


# In[4]:


lytic_genes = df.loc[df['concat_annotations'].str.contains(lytic_keywords)]


# In[5]:


lytic_genes['gene_position'] = lytic_genes['gene_position'].astype(int)


# In[6]:


# Add gene key so we can filter the htseq table
gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt', sep='\t')
gff_df = gff_df[['vOTU', 'gene']]
gff_df['gene_position'] = gff_df['gene'].apply(lambda x: int(x.split('_')[1]))
gff_df.head()


# In[7]:


lytic_genes = lytic_genes.merge(gff_df, on=['vOTU', 'gene_position'], how='left')


# In[8]:


# number of vOTUs with a lytic gene
lytic_genes['vOTU'].nunique(), lytic_genes['vOTU'].nunique()/df.vOTU.nunique()


# 77% of vOTUs identified have at least one of these genes

# In[9]:


# get the gene IDs of each gene with a structural phage or lytic gene annotation
lytic_gene_list = list(lytic_genes['gene'].unique())


# ## 2. Load htseq count table and set all read counts <5 to zero.

# In[10]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import observed_otus
import seaborn as sns

# Read data
htseq_counts_path = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED.tsv'

data = pd.read_csv(htseq_counts_path, sep='\t', header=None)

# Assign column names. These are in the correct order from htseq counts but no columns are assigned.
columns = ["gene",
  "STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004",
  "STM_0716_E_M_E025", "STM_0716_E_M_E027", "STM_0716_E_M_E029",
  "STM_0716_E_M_E030", "STM_0716_E_M_E031", "STM_0716_E_M_E033",
  "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
  "STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E054",
  "STM_0716_E_M_E055", "STM_0716_E_M_E056", "STM_0716_E_M_E058",
  "STM_0716_E_M_E059", "STM_0716_E_M_E060", "STM_0716_E_M_E062",
  "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E066",
  "STM_0716_E_M_E067", "STM_0716_E_M_E068", "STM_0716_E_M_E070",
  "STM_0716_E_M_E071", "STM_0716_E_M_E072", "STM_0716_E_M_E121",
  "STM_0716_E_M_E122", "STM_0716_E_M_E123", "STM_0716_E_M_E125",
  "STM_0716_E_M_E126", "STM_0716_E_M_E127", "STM_0716_E_M_E129",
  "STM_0716_E_M_E130", "STM_0716_E_M_E131"]

data.columns = columns

# Drop condensed tannin sample columns
data = data.drop(columns=['STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                       'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                       'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                       'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127'])


# In[11]:


data


# In[12]:


# graph gene mapping stats

# get lib size
lib = pd.DataFrame(data.iloc[:, 1:27].sum()).reset_index()
lib.columns = ['sample', 'total']

# Get no_feature and ambig counts (last two rows)
lib['no_feature'] = data.iloc[-5, 1:27].values
lib['ambig'] = data.iloc[-4, 1:27].values

# Calculate mapped reads
lib['mapped'] = lib['total'] - lib['no_feature'] - lib['ambig']

# Reshape data for plotting
lib_melted = lib.melt(id_vars=['sample', 'total'], 
                      value_vars=['no_feature', 'ambig', 'mapped'],
                      var_name='status', value_name='counts')
lib_melted


# In[13]:


lib_melted.loc[lib_melted['sample'] == 'STM_0716_E_M_E035']


# In[14]:


# Alternative using pandas plotting (simpler but less customizable)
lib[['mapped', 'ambig', 'no_feature']].plot.barh(
    stacked=True, 
    figsize=(12, 8),
    color=['#1f77b4', '#ff7f0e', '#2ca02c']
)

plt.xlabel('Counts')
plt.ylabel('Sample')
plt.title('Mapping Statistics - Stacked')
plt.legend(title='Status', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# In[15]:


# Drop last five metadata rows
data = data.iloc[:-5, :]

# Add vOTU key merging on gene ID
data = data.merge(gff_df, on='gene', how='left')
data['is_lytic'] = data['gene'].isin(lytic_gene_list)

data.tail()


# In[16]:


import numpy as np
import matplotlib.pyplot as plt

def plot_dataframe_histogram(df):
    """
    Plot histogram of non-zero values and print proportion of zeros.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to analyze
    """
    # Flatten all values
    all_values = df.values.flatten()
    
    # Calculate proportion of zeros
    zero_count = np.sum(all_values == 0)
    total_count = len(all_values)
    zero_proportion = zero_count / total_count
    
    print(f"Proportion of gene-sample pairs with a zero value: {zero_proportion:.3f} ({zero_count}/{total_count})\n")
    
    # Get non-zero values
    nonzero_values = all_values[all_values != 0]

    percentiles = np.arange(0, 101, 10)  # 0, 10, 20, ..., 100
    values = np.percentile(nonzero_values, percentiles)
    for p, v in zip(percentiles, values):
        print(f"{p:3d}th percentile: {v:.6f}")
    
    if len(nonzero_values) == 0:
        print("No non-zero values to plot")
        return
    
    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.hist(nonzero_values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.yscale('log')
    plt.title('Histogram of number of reads mapping to genes')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

# Example usage:
plot_dataframe_histogram(data.iloc[:, 1:27])


# In[17]:


# Make a copy to mask such that counts < 5 are set to zero
data_5x = data.copy()

# Select sample columns
numeric_cols = data_5x.columns[1:27]

# Step 3: Set values < 5 to 0 in those columns
data_5x[numeric_cols] = data_5x[numeric_cols].apply(lambda x: x.mask(x < 5, 0))


# In[18]:


data_5x


# In[19]:


data_5x_melted = data_5x.melt(id_vars=['vOTU', 'gene', 'gene_position', 'is_lytic'], var_name='Sample', value_name='readcount')
data_5x_melted.head()


# In[20]:


data_5x_melted['is_lytic_and_readcount_5x'] = (data_5x_melted['is_lytic'] == True) & (data_5x_melted['readcount'] >= 5)


# In[21]:


vOTU_lytic_in_sample = data_5x_melted.groupby(['vOTU', 'Sample'])['is_lytic_and_readcount_5x'].any().reset_index()
vOTU_lytic_in_sample.head()


# In[22]:


# Add conditional to figure out if vOTUs are covered at least 10% on the reverse strand
prop_covered_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/prop_covered_metaT/all_proportions_reverse.txt', sep='\t', header=None)
prop_covered_df.columns = ['vOTU', 'Sample', 'prop_covered']
prop_covered_df['is_covered_10perc'] = prop_covered_df['prop_covered'] >= 0.1

# merge and fill NA with 0 because NA in prop covered means no reads mapped in that sample to that contig.
vOTU_lytic_in_sample_prop_covered = vOTU_lytic_in_sample.merge(prop_covered_df[['vOTU', 'Sample', 'prop_covered']], on=['vOTU', 'Sample'], how='left').fillna(0)

# create new conditional. We will use this table later to filter GeTMM relative abundance values.
vOTU_lytic_in_sample_prop_covered['is_active'] = (vOTU_lytic_in_sample_prop_covered['is_lytic_and_readcount_5x'] == True) & (vOTU_lytic_in_sample_prop_covered['prop_covered'] >= 0.1)
vOTU_lytic_in_sample_prop_covered.head()


# ## Compute GeTMM relative abundance of each vOTU

# In[23]:


# get read counts in reads per 1000 (rpk)
len_df = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt", sep="\t")
len_df['length'] = len_df['length'].astype(int)
data_5x_len = data_5x.merge(len_df, on=['vOTU','gene'], how='left')

# Make sure gene is index
data_5x_len = data_5x_len.set_index('gene')

# Reindex lengths so they match
len_kb = (data_5x_len['length'] / 1000).reindex(data_5x_len.index)

# Divide (row-wise, with proper alignment)
data_5x_rpk = data_5x_len.iloc[:, :26].divide(len_kb, axis=0)

data_5x_rpk.head()


# In[24]:


data_5x_rpk.to_csv('data_5x_rpk.csv')


# In[25]:


# Compute getmms in EdgeR by feeding data_5x_rpk.csv to getmm.R
# I tried using conorm, a package in python for getmm, but I do not trust the TMM algorithm as it produces different results to the EdgeR implementation.
'''
library(edgeR)
rpk = read.delim('/users/PAS1573/riddell26/data_5x_rpk.csv', sep=',')
rownames(rpk) <- rpk$gene
rpk$gene <- NULL
group <- c(rep("A",ncol(rpk)))

rpk.norm <- DGEList(counts=rpk,group=group)
rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df)
# 0 1142551
'''


# In[26]:


getmms_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs.tsv', sep='\t')
getmms_df = getmms_df.reset_index(names='gene').merge(data_5x[['vOTU', 'gene']], how='left')
getmms_df.head()


# In[27]:


# Get the mean getmm value per vOTU in each sample
getmms_vOTU = (
    getmms_df
    .set_index('gene')
    .groupby('vOTU')
    .mean()
    .reset_index()
    .melt(id_vars='vOTU', var_name='Sample', value_name='getmm_abundance')
)
getmms_vOTU


# In[28]:


# merge with vOTU_lytic_in_sample_prop_covered[['vOTU', 'Sample' 'is_active']] and set to zero if is_active is False
getmms_vOTU_corrected = getmms_vOTU.merge(vOTU_lytic_in_sample_prop_covered[['vOTU', 'Sample', 'is_active']], on=['vOTU', 'Sample'])
getmms_vOTU_corrected['getmm_abundance_corrected'] = getmms_vOTU_corrected['getmm_abundance'].where(getmms_vOTU_corrected['is_active'], 0)
getmms_vOTU_corrected


# In[29]:


# write to csv
getmms_vOTU_corrected[['vOTU', 'Sample', 'getmm_abundance_corrected']].to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_melted.tsv', sep='\t')

# write unmelted as well
getmms_vOTU_corrected_wide = getmms_vOTU_corrected.pivot(
    index="vOTU", 
    columns="Sample", 
    values="getmm_abundance_corrected"
)

# drop if all zeroes
getmms_vOTU_corrected_wide = getmms_vOTU_corrected_wide.loc[getmms_vOTU_corrected_wide.sum(axis=1) > 0]
print('num vOTUs active in at least one sample: ' + str(len(getmms_vOTU_corrected_wide)))

getmms_vOTU_corrected_wide.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')


# In[32]:


# 1. Convert active_df to binary (1 = active, 0 = inactive)
binary_df = (getmms_vOTU_corrected_wide > 0).astype(int)

# 2. Count number of samples each vOTU was active in (row-wise sum)
votu_activity_counts = binary_df.sum(axis=1)

# 3. Count number of active vOTUs in each sample (column-wise sum)
sample_activity_counts = binary_df.sum(axis=0)

# Create subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# 4. Left plot: how many samples each vOTU was active in
sns.violinplot(y=votu_activity_counts, ax=ax1)
ax1.axhline(y=np.median(votu_activity_counts), color='red', linestyle='--', linewidth=2, 
            label=f'Median: {np.median(votu_activity_counts):.1f}')
ax1.set_title("Number of Samples Each vOTU Was Active In")
ax1.set_ylabel("Number of Samples")
ax1.set_xlabel("vOTUs")
ax1.legend()

# 5. Right plot: how many vOTUs were active in each sample
sns.violinplot(y=sample_activity_counts, ax=ax2)
ax2.axhline(y=np.median(sample_activity_counts), color='red', linestyle='--', linewidth=2, 
            label=f'Median: {np.median(sample_activity_counts):.1f}')
ax2.set_title("Number of Active vOTUs per Sample")
ax2.set_ylabel("Number of Active vOTUs")
ax2.set_xlabel("Samples")
ax2.legend()

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1E_vOTUs_vs_samples.pdf', dpi=300)
plt.show()


# In[ ]:




