#!/usr/bin/env python
# coding: utf-8

# # This notebook plots alpha diversity statistics based on vOTU GeTMM relative abundance computed from metatranscriptomes

# In[1]:


import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import regex as re
import textwrap
import time
from matplotlib_venn import venn3
import skbio
import matplotlib as mpl


# In[2]:


from print_versions import print_versions
print_versions(globals())


# In[3]:


# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


# In[4]:


# df generated from 06-identify-active-vOTU-and-calc-getmm-abundance.py
df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
df.head()


# In[5]:


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


# In[6]:


replicate_frame = replicate_frame.loc[replicate_frame['treatment'] != 'CT']


# In[7]:


replicate_frame.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv', index=False)


# # Richness

# In[8]:


df = df.set_index('vOTU')


# In[9]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import observed_otus
import seaborn as sns

# Read data
df = np.ceil(df)
df.head()


# In[10]:


# Read metadata
metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")

# Simplify treatment labels
metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
metadata = metadata.loc[metadata['treatment'] != 'T']

# Create treatment_timepoint column
treatment_timepoint = metadata['treatment'] + "_" + metadata['timepoint'].astype(str)
df.columns = treatment_timepoint.values

# Make unique column names
def make_unique(columns):
    seen = {}
    new_cols = []
    for col in columns:
        if col not in seen:
            seen[col] = 0
            new_cols.append(col)
        else:
            seen[col] += 1
            new_cols.append(f"{col}.{seen[col]}")
    return new_cols

df.columns = make_unique(df.columns)

df


# In[11]:


# Transpose and convert to numeric matrix
data_mtx = df.T
data_mtx = data_mtx.apply(pd.to_numeric, errors='coerce')
data_mtx


# In[12]:


# Subsample and compute rarefied richness
S = data_mtx.astype(int).apply(observed_otus, axis=1)


# In[13]:


S


# In[14]:


metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('day')[1]))


# In[15]:


metadata


# In[16]:


metadata['richness'] = S.values


# In[17]:


metadata


# In[18]:


# Plot individual points and a mean line per treatment
plt.figure(figsize=(6, 6))

# Points (individual replicates)
sns.scatterplot(
    data=metadata,
    x='timepoint',
    y='richness',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    alpha=1
)

# Mean line across timepoints per treatment
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='richness',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=4,
    legend=False  # Prevent duplicate legends
)

# plt.title("Richness over Time by Treatment")
plt.xlabel("Day", fontsize=14)
plt.ylabel("Richness", fontsize=14)
treatment_colors = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

plt.xticks([0, 7, 14, 21, 35], fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, alpha=0.3)

# Remove legend
plt.legend().remove()

# Remove plot borders (spines)
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_richness.pdf', dpi=300)
plt.show()



# # Shannon Diversity

# In[19]:


from skbio.diversity.alpha import shannon
data_mtx.apply(shannon, axis=1)


# In[20]:


metadata['shannon_diversity'] = data_mtx.apply(shannon, axis=1).values


# In[21]:


# Plot individual points and a mean line per treatment
plt.figure(figsize=(6, 6))

# Points (individual replicates)
sns.scatterplot(
    data=metadata,
    x='timepoint',
    y='shannon_diversity',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    alpha=1
)

# Mean line across timepoints per treatment
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='shannon_diversity',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=4,
    legend=False  # Prevent duplicate legends
)

# plt.title("Shannon Diversity")
plt.xlabel("Day", fontsize=14)
plt.ylabel("Shannon Index", fontsize=14)
treatment_colors = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

plt.xticks([0, 7, 14, 21, 35], fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, alpha=0.3)

# Remove legend
plt.legend().remove()

# Remove plot borders (spines)
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_shannon_diversity.pdf', dpi=300)
plt.show()



# # Inverse Simpson

# In[22]:


from skbio.diversity.alpha import inv_simpson
metadata['inverse_simpson'] = data_mtx.apply(lambda x: inv_simpson(x), axis=1).values


# In[23]:


metadata


# In[24]:


# Plot individual points and a mean line per treatment
plt.figure(figsize=(6, 6))

# Points (individual replicates)
sns.scatterplot(
    data=metadata,
    x='timepoint',
    y='inverse_simpson',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    alpha=1
)

# Mean line across timepoints per treatment
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='inverse_simpson',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=4,
    legend=False  # Prevent duplicate legends
)

# plt.title("Simpson Diversity Index (1-D)")
plt.xlabel("Day", fontsize=14)
plt.ylabel("Inverse Simpson's Index (1/D)", fontsize=14)
treatment_colors = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

plt.xticks([0, 7, 14, 21, 35], fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, alpha=0.3)

# Remove legend
plt.legend().remove()

# Remove plot borders (spines)
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_inverse_simpson_diversity.pdf', dpi=300)
plt.show()


# # Proportion of reads mapping to active vOTUs

# In[25]:


gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt', sep='\t')
gff_df = gff_df[['vOTU', 'gene']]
gff_df.head()

# Read the TSV file without header
counts = pd.read_csv(
    "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED.tsv",
    sep="\t",
    header=None
)

# Assign column names
counts.columns = [
    "gene",
    "STM_0716_E_M_E002",
    "STM_0716_E_M_E003",
    "STM_0716_E_M_E004",
    "STM_0716_E_M_E025",
    "STM_0716_E_M_E027",
    "STM_0716_E_M_E029",
    "STM_0716_E_M_E030",
    "STM_0716_E_M_E031",
    "STM_0716_E_M_E033",
    "STM_0716_E_M_E034",
    "STM_0716_E_M_E035",
    "STM_0716_E_M_E050",
    "STM_0716_E_M_E051",
    "STM_0716_E_M_E052",
    "STM_0716_E_M_E054",
    "STM_0716_E_M_E055",
    "STM_0716_E_M_E056",
    "STM_0716_E_M_E058",
    "STM_0716_E_M_E059",
    "STM_0716_E_M_E060",
    "STM_0716_E_M_E062",
    "STM_0716_E_M_E063",
    "STM_0716_E_M_E064",
    "STM_0716_E_M_E066",
    "STM_0716_E_M_E067",
    "STM_0716_E_M_E068",
    "STM_0716_E_M_E070",
    "STM_0716_E_M_E071",
    "STM_0716_E_M_E072",
    "STM_0716_E_M_E121",
    "STM_0716_E_M_E122",
    "STM_0716_E_M_E123",
    "STM_0716_E_M_E125",
    "STM_0716_E_M_E126",
    "STM_0716_E_M_E127",
    "STM_0716_E_M_E129",
    "STM_0716_E_M_E130",
    "STM_0716_E_M_E131"
]
counts.head()


# In[26]:


counts = counts.merge(gff_df, on='gene', how='left').dropna()
counts


# In[27]:


counts_long = counts.melt(id_vars=['gene', 'vOTU'], var_name='Sample', value_name='num_reads_mapped')
counts_long


# In[28]:


counts_long = counts_long.groupby(['vOTU', 'Sample']).agg({'num_reads_mapped': 'sum'}).reset_index()


# In[29]:


counts_long = counts_long.merge(replicate_frame, on='Sample', how='left')


# In[30]:


counts_long = counts_long.loc[counts_long['treatment'] != 'CT']


# In[31]:


counts_long


# In[32]:


# keep only active vOTUs
counts_long = counts_long.loc[counts_long['vOTU'].isin(list(df.index))]
counts_long


# In[33]:


reads_mapped_per_sample = counts_long.groupby(['Sample', 'treatment', 'day', 'replicate']).agg({'num_reads_mapped': 'sum'}).reset_index()
reads_mapped_per_sample


# In[34]:


metadata.columns


# In[35]:


reads_mapped_per_sample = reads_mapped_per_sample.merge(metadata[['Sample', 'total reads (R1+R2)']], on='Sample', how='left')


# In[36]:


reads_mapped_per_sample['prop_mapped'] = reads_mapped_per_sample['num_reads_mapped'] / reads_mapped_per_sample['total reads (R1+R2)'] * 100
reads_mapped_per_sample


# In[37]:


reads_mapped_per_sample['day'] = reads_mapped_per_sample['day'].astype(int)


# In[38]:


reads_mapped_per_sample.loc[reads_mapped_per_sample['day'] == 14]


# In[39]:


metadata = metadata.merge(reads_mapped_per_sample[['Sample', 'prop_mapped']], on='Sample', how='left')
metadata


# In[43]:


import matplotlib.pyplot as plt
import seaborn as sns

# Base plot
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot shannon_diversity on the primary y-axis
# sns.scatterplot(
#     data=metadata,
#     x='timepoint',
#     y='shannon_diversity',
#     hue='treatment',
#     palette={"C": "#FC9B2D", "U": "#7ACBC3"},
#     alpha=1,
#     ax=ax1
# )
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='shannon_diversity',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=4,
    legend=False,
    ax=ax1
)

ax1.set_xlabel("Day", fontsize=14)
ax1.set_ylabel("vOTU Shannon Index (metaT)", fontsize=14)
ax1.set_xticks([0, 7, 14, 21, 35])
ax1.tick_params(axis='both', labelsize=12)
ax1.grid(True, alpha=0.3)

# Remove plot borders
for spine in ax1.spines.values():
    spine.set_visible(False)

# Create second y-axis
ax2 = ax1.twinx()
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='prop_mapped',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    linestyle='--',
    errorbar='ci',
    linewidth=4,
    legend=False,
    ax=ax2
)
ax2.set_ylabel("% of total RNA mapped to active vOTUs", fontsize=14)
ax2.tick_params(axis='y', labelsize=12)

# Final adjustments
fig.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-A_shannon_diversity_with_prop_mapped.pdf', dpi=300)
plt.show()


# In[41]:


import matplotlib.pyplot as plt
import seaborn as sns

# Base plot
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot shannon_diversity on the primary y-axis
# sns.scatterplot(
#     data=metadata,
#     x='timepoint',
#     y='shannon_diversity',
#     hue='treatment',
#     palette={"C": "#FC9B2D", "U": "#7ACBC3"},
#     alpha=1,
#     ax=ax1
# )
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='shannon_diversity',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=4,
    legend=False,
    ax=ax1
)

ax1.set_xlabel("Day", fontsize=14)
ax1.set_ylabel("vOTU Shannon Index (metaT)", fontsize=18)
ax1.set_xticks([0, 7, 14, 21, 35])
ax1.tick_params(axis='both', labelsize=18)
ax1.grid(True, alpha=0.3)

# Remove plot borders
for spine in ax1.spines.values():
    spine.set_visible(False)

# Create second y-axis
ax2 = ax1.twinx()
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='prop_mapped',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    linestyle='--',
    errorbar='ci',
    linewidth=4,
    legend=False,
    ax=ax2
)
ax2.set_ylabel("% of total RNA mapped\nto active vOTUs", fontsize=18)
ax2.tick_params(axis='y', labelsize=18)

# Final adjustments
fig.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-A_shannon_diversity_with_prop_mapped_for_presentation.pdf', dpi=300)
plt.show()


# In[42]:


get_ipython().system('jupyter nbconvert --to script 003-alpha-diversity.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/scripts/07-alpha_diversity')


# In[ ]:




