#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


# from predicted_hosts_EDA.ipynb
iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/results/iphop_custom/Host_prediction_to_genome_m90.csv')
iphop_df.head()


# In[3]:


# load in updated host classification from GTDB-classify
arch = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/data/gtdb_classify_wf_for_virmatcher/gtdbtk.ar122.summary.tsv', sep='\t')
arch['Host genome'] = arch['user_genome']
arch['Host taxonomy'] = arch['classification']
arch = arch[['Host genome', 'Host taxonomy']]

bac = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/data/gtdb_classify_wf_for_virmatcher/gtdbtk.bac120.summary.tsv', sep='\t')
bac['Host genome'] = bac['user_genome']
bac['Host taxonomy'] = bac['classification']
bac = bac[['Host genome', 'Host taxonomy']]

host_taxonomy = pd.concat([arch, bac])


# In[4]:


# replace host taxonomy for hosts from the custom iphop database

# Merge to get updated taxonomy only for matching 'Add_' host genomes
iphop_df = iphop_df.merge(host_taxonomy, on="Host genome", how="left", suffixes=("", "_new"))

# Replace 'Host taxonomy' where 'Host genome' starts with 'Add_'
iphop_df["Host taxonomy"] = iphop_df.apply(
    lambda row: row["Host taxonomy_new"] if row["Host genome"].startswith("Add_") and pd.notna(row["Host taxonomy_new"]) else row["Host taxonomy"],
    axis=1
)

# Drop the extra merged column
iphop_df.drop(columns=["Host taxonomy_new"], inplace=True)


# In[17]:


# How many vOTUs had an iPHoP prediction?=
print('number of vOTUs from reference dataset with a host prediction: ' + str(iphop_df.Virus.nunique()))


# In[16]:


# How many transciptionally active vOTUs had an iPHoP prediction?
with open('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/results/metaT/active_1_gene_per_10kb_vOTUs.txt', 'r') as f:
    active_vOTUs = [i.strip() for i in f.readlines()]

num_active_predictions = iphop_df.loc[iphop_df['Virus'].isin(active_vOTUs)].Virus.nunique()
print('number of active vOTUs with a host prediction: ' + str(num_active_predictions))


# In[23]:


# number of iphop predictions with a Stordalen host prediction
print('number of vOTUs from reference dataset with a host prediction: ' + str(iphop_df.loc[iphop_df['Host genome'].str.startswith('Add_')].Virus.nunique()))


# In[25]:


# number of iphop predictions with a Stordalen host prediction
print('number of active vOTUs from reference dataset with a host prediction: ' + str(iphop_df
                                                                                     .loc[
                                                                                     (iphop_df['Host genome'].str.startswith('Add_'))
                                                                                     & (iphop_df['Virus'].isin(active_vOTUs)) 
                                                                                      ].Virus.nunique()))


# In[34]:


# number of active vOTUs with an active host
MAG_abundance = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/meancov/metaG/MAGs/metaG_MAG_mean_97_75_75.tsv', sep='\t')
present_MAGs = MAG_abundance.loc[MAG_abundance.iloc[:,1:].sum(axis=1) > 0]


# In[38]:


present_MAGs_list = ['Add_' + i for i in present_MAGs['Genome']]


# In[39]:


active_vOTUs_present_MAGs = (iphop_df
     .loc[
     (iphop_df['Host genome'].isin(present_MAGs_list))
     & (iphop_df['Virus'].isin(active_vOTUs)) 
      ]
)


# In[41]:


len(active_vOTUs_present_MAGs)


# In[42]:


# build virus-host matrix with relevant metadata


# In[48]:


# add provirus
vOTU_metadata = pd.read_csv('/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
is_provirus = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()


# In[50]:


iphop_df = iphop_df.rename(columns={'Virus': 'vOTU'})


# In[51]:


iphop_df = iphop_df.merge(is_provirus, on='vOTU', how='left')



# In[54]:


iphop_df[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = iphop_df['Host taxonomy'].str.split(';', expand=True)


# In[55]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

iphop_df['highest_host_tax_rank'] = iphop_df.apply(get_highest_resolution, axis=1)


# In[63]:


iphop_df['vOTU_is_active'] = iphop_df['vOTU'].isin(active_vOTUs)
iphop_df['MAG_is_present'] = iphop_df['Host genome'].isin(present_MAGs_list)
iphop_df['MAG_is_stordalen'] = iphop_df['Host genome'].str.startswith('Add_')


# In[65]:


iphop_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_for_cytoscape.tsv', sep='\t', index=False)


# In[66]:


active_vOTUs_present_MAGs = (iphop_df
     .loc[
            (iphop_df['vOTU_is_active'] == True) 
          & (iphop_df['MAG_is_present'] == True)
         ]
)


# In[68]:


active_vOTUs_present_MAGs.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/active_vOTUs_present_MAGs_for_cytoscape.tsv', sep='\t', index=False)

