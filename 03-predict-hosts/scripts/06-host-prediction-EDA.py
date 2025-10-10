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
arch = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/data/gtdb_classify_wf/gtdbtk.ar122.summary.tsv', sep='\t')
arch['Host genome'] = arch['user_genome']
arch['Host taxonomy'] = arch['classification']
arch = arch[['Host genome', 'Host taxonomy']]

bac = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/data/gtdb_classify_wf/gtdbtk.bac120.summary.tsv', sep='\t')
bac['Host genome'] = bac['user_genome']
bac['Host taxonomy'] = bac['classification']
bac = bac[['Host genome', 'Host taxonomy']]

host_taxonomy = pd.concat([arch, bac])
host_taxonomy.head()


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
iphop_df.head()


# In[5]:


iphop_df['Virus'].nunique()


# # Stitch: extract ALL vOTUs, active or inactive, that were predicted to infect JAGFXR01s. Need it for figure 3 comparative genomics

# In[6]:


jag_viruses = list(iphop_df.loc[iphop_df['Host taxonomy'].str.contains('JAGFXR01')]['Virus'].unique())


# In[7]:


with open('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/JAGFXR01_vOTUs/JAGFXR01_vOTUs.txt', 'w') as f:
    for i in jag_viruses:
        f.write(i)
        f.write('\n')


# In[8]:


iphop_df.loc[iphop_df['Virus'] == 'contig_1090030']


# # Stitch #2, investigate if top vOTUs had ANY host prediction

# In[9]:


# from 004-vOTU-relative-abundance-over-time-getmm-by-vOTU, df.loc[(df['sum_getmm'] > 11749.564422) & (df['highest_host_tax_rank'] == 'no_assignment')]['vOTU'].unique()
top_vOTUs = [
    'contig_1022629',
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


# In[10]:


iphop_df.loc[iphop_df['Virus'].isin(top_vOTUs)]


# # Continue filtering for only active vOTUs

# In[11]:


# How many vOTUs had an iPHoP prediction?
print('number of vOTUs from reference dataset with a host prediction: ' + str(iphop_df.Virus.nunique()))


# In[12]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
active_vOTUs = list(df['vOTU'].unique())


# In[13]:


num_active_predictions = iphop_df.loc[iphop_df['Virus'].isin(active_vOTUs)].Virus.nunique()
print('number of active vOTUs with a host prediction: ' + str(num_active_predictions))


# In[14]:


# number of iphop predictions with a Stordalen host prediction
print('number of vOTUs from reference dataset with a bioreactor host prediction: ' + str(iphop_df.loc[iphop_df['Host genome'].str.startswith('Add_')].Virus.nunique()))


# In[15]:


# number of iphop predictions with a Stordalen host prediction
print('number of active vOTUs from reference dataset with a bioreactor host prediction: ' + str(iphop_df
                                                                                     .loc[
                                                                                     (iphop_df['Host genome'].str.startswith('Add_'))
                                                                                     & (iphop_df['Virus'].isin(active_vOTUs)) 
                                                                                      ].Virus.nunique()))


# In[16]:


# number of active vOTUs with an active host
active_MAGs = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/geTMM_table.csv')
len(active_MAGs)


# In[17]:


active_MAGs = active_MAGs.groupby('fasta').sum()


# In[18]:


active_MAGs = active_MAGs.reset_index().rename(columns={'fasta': 'MAG'})
active_MAGs.head()


# In[19]:


len(active_MAGs)


# In[20]:


active_MAGs


# In[21]:


with open('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_list_no_prefix.txt', 'w') as f:
    for i in list(active_MAGs['MAG']):
        f.write(i)
        f.write('\n')


# In[22]:


active_MAGs_list = set(['Add_' + i for i in active_MAGs['MAG']])


# In[23]:


len(active_MAGs_list)


# In[24]:


with open('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_list.txt', 'w') as f:
    for i in list(active_MAGs_list):
        f.write(i)
        f.write('\n')


# In[25]:


active_vOTUs_active_MAGs = (iphop_df
     .loc[
     (iphop_df['Host genome'].isin(active_MAGs_list))
     & (iphop_df['Virus'].isin(active_vOTUs)) 
      ]
)


# In[26]:


len(active_vOTUs_active_MAGs)


# In[27]:


active_vOTUs_active_MAGs.Virus.nunique()


# In[28]:


active_vOTUs_active_MAGs['Host genome'].nunique()


# # wrangle provirus metadata to add to iPHoP results for cytoscape visualization

# In[30]:


vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
is_provirus = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()


# In[31]:


vOTU_metadata.columns


# In[32]:


# export provirus and their coordinates for induced prophage analysis
prophage_df = vOTU_metadata.loc[vOTU_metadata['is_provirus_aggregated'] == True][['vOTU', 'is_provirus_aggregated', 'coordinates']].drop_duplicates()
prophage_df


# In[33]:


genomad_prophage = prophage_df.loc[prophage_df['vOTU'].str.contains('provirus')]
genomad_prophage.head()


# In[34]:


genomad_prophage['scaffold'] = genomad_prophage['vOTU'].apply(lambda x: x.split('|')[0])


# In[35]:


genomad_prophage[['start', 'stop']] = genomad_prophage['coordinates'].str.split('-', expand=True)


# In[36]:


genomad_prophage


# In[37]:


genomad_prophage = genomad_prophage[['scaffold', 'vOTU', 'start', 'stop']]
genomad_prophage.columns = ['scaffold', 'fragment', 'start', 'stop']


# In[38]:


# make sure contig # is before >contig_432973 since these are what are on the MAGs.
genomad_prophage = genomad_prophage.loc[genomad_prophage['scaffold'].apply(lambda x: int(x.split('_', 1)[1]) <= 432973)]


# In[39]:


genomad_prophage.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/09-predict-prophage-activity/data/prophage_coordinates_for_propagate.tsv', sep='\t', index=False)


# # Merge with iPHoP

# In[40]:


iphop_df = iphop_df.rename(columns={'Virus': 'vOTU'})


# In[41]:


iphop_df = iphop_df.merge(is_provirus, on='vOTU', how='left')


# In[42]:


iphop_df


# In[43]:


iphop_df[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = iphop_df['Host taxonomy'].str.split(';', expand=True)


# In[44]:


iphop_df['G'].nunique()


# In[45]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

iphop_df['highest_host_tax_rank'] = iphop_df.apply(get_highest_resolution, axis=1)


# In[46]:


iphop_df['vOTU_is_active'] = iphop_df['vOTU'].isin(active_vOTUs)
iphop_df['MAG_is_active'] = iphop_df['Host genome'].isin(active_MAGs_list)
iphop_df['MAG_is_stordalen'] = iphop_df['Host genome'].str.startswith('Add_')


# In[47]:


iphop_df.head()


# In[48]:


iphop_df.loc[iphop_df['vOTU'].isin(top_vOTUs)].drop_duplicates(subset=['vOTU'])


# In[49]:


iphop_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t', index=False)


# # Visualize host prediction stats

# In[51]:


iphop_host_stats = iphop_df.drop_duplicates(['Host genome']).groupby(['highest_host_tax_rank']).agg({'Host genome': 'nunique', 'MAG_is_active': 'sum'}).reset_index().sort_values(by='MAG_is_active', ascending=False)
iphop_host_stats


# In[52]:


iphop_virus_stats = iphop_df.drop_duplicates(['vOTU', 'highest_host_tax_rank']).groupby(['highest_host_tax_rank']).agg({'vOTU': 'nunique', 'vOTU_is_active': 'sum'}).reset_index().sort_values(by='vOTU', ascending=False)


# In[53]:


iphop_virus_stats


# In[54]:


combined_stats = iphop_host_stats.merge(iphop_virus_stats, on='highest_host_tax_rank', how='outer')
combined_stats


# In[55]:


combined_stats.loc[combined_stats['highest_host_tax_rank'] == 'g__Terracidiphilus']


# In[56]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Assuming your data is in a DataFrame called 'df'

# Sort by vOTU count in descending order
df_sorted = combined_stats.sort_values('vOTU_is_active', ascending=False)

# Get top 20 categories by vOTU count
top_20_votu = df_sorted.head(20)

# Calculate "Other" category for vOTU and vOTU_is_active
other_votu = df_sorted.iloc[20:]['vOTU'].sum()
other_votu_active = df_sorted.iloc[20:]['vOTU_is_active'].sum()

# Create DataFrame for vOTU chart with scaled "Other"
votu_chart_data = top_20_votu[['highest_host_tax_rank', 'vOTU', 'vOTU_is_active']].copy()
votu_scale_factor = 0.3  # Scale down the "Other" bar for vOTU
votu_chart_data = pd.concat([votu_chart_data, 
                           pd.DataFrame([['Other', other_votu * votu_scale_factor, other_votu_active * votu_scale_factor]], 
                                       columns=['highest_host_tax_rank', 'vOTU', 'vOTU_is_active'])])

# Store actual values for labeling
votu_actual_values = top_20_votu[['vOTU', 'vOTU_is_active']].copy()
votu_actual_values = pd.concat([votu_actual_values, 
                              pd.DataFrame([[other_votu, other_votu_active]], 
                                          columns=['vOTU', 'vOTU_is_active'])])

# Sort by MAG count in descending order for MAG chart
df_sorted_mag = combined_stats.sort_values('MAG_is_active', ascending=False)
top_20_mag = df_sorted_mag.head(20)

# Calculate "Other" category for MAG and MAG_is_active
other_mag = df_sorted_mag.iloc[20:]['Host genome'].sum()
other_mag_active = df_sorted_mag.iloc[20:]['MAG_is_active'].sum()

# Create DataFrame for MAG chart with scaled "Other"
mag_scale_factor = 0.3  # Scale down the "Other" bar for MAG
mag_chart_data = top_20_mag[['highest_host_tax_rank', 'Host genome', 'MAG_is_active']].copy()
mag_chart_data = pd.concat([mag_chart_data, 
                          pd.DataFrame([['Other', other_mag * mag_scale_factor, other_mag_active * mag_scale_factor]], 
                                      columns=['highest_host_tax_rank', 'Host genome', 'MAG_is_active'])])

# Store actual values for labeling
mag_actual_values = top_20_mag[['Host genome', 'MAG_is_active']].copy()
mag_actual_values = pd.concat([mag_actual_values, 
                             pd.DataFrame([[other_mag, other_mag_active]], 
                                         columns=['Host genome', 'MAG_is_active'])])

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Plot vOTU chart with stacked bars
x_pos = range(len(votu_chart_data))

# Plot total vOTU (white bars)
ax1.bar(x_pos, votu_chart_data['vOTU'], color='white', edgecolor='black', linewidth=0.5, label='Total vOTU')

# Plot active vOTU on top (black bars)
ax1.bar(x_pos, votu_chart_data['vOTU_is_active'], color='black', label='Active vOTU')

# ax1.set_title('vOTU Distribution (Top 20 + Other)', fontsize=14, fontweight='bold')
ax1.set_xlabel('Host Taxonomy Rank', fontsize=12)
ax1.set_ylabel('vOTU Count', fontsize=12)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(votu_chart_data['highest_host_tax_rank'], rotation=45, ha='right')
ax1.legend()

# Add value labels on vOTU bars (using actual values, not scaled)
for i, (scaled_val, actual_val) in enumerate(zip(votu_chart_data['vOTU'], votu_actual_values['vOTU'])):
    if votu_chart_data['highest_host_tax_rank'].iloc[i] == 'Other':
        # Position label above the scaled bar for "Other"
        ax1.text(i, scaled_val + max(votu_chart_data['vOTU'][:-1]) * 0.05, str(int(actual_val)), 
                ha='center', va='bottom', fontsize=9, fontweight='bold', color='red')
    else:
        ax1.text(i, scaled_val + 0.1, str(int(actual_val)), ha='center', va='bottom', fontsize=9)

# Add value labels for active vOTU (using actual values)
for i, (scaled_total, scaled_active, actual_active) in enumerate(zip(
    votu_chart_data['vOTU'], votu_chart_data['vOTU_is_active'], votu_actual_values['vOTU_is_active'])):
    
    if actual_active > 0:
        if votu_chart_data['highest_host_tax_rank'].iloc[i] == 'Other':
            ax1.text(i, scaled_active/2, str(int(actual_active)), ha='center', va='center', 
                    color='pink', fontsize=9, fontweight='bold')
        else:
            ax1.text(i, scaled_active/2, str(int(actual_active)), ha='center', va='center', 
                    color='pink', fontsize=9, fontweight='bold')

# Add break indicator for "Other" in vOTU chart
other_idx_votu = len(votu_chart_data) - 1
ax1.plot([other_idx_votu-0.4, other_idx_votu+0.4], 
        [votu_chart_data['vOTU'].iloc[other_idx_votu] + max(votu_chart_data['vOTU'][:-1]) * 0.02, 
         votu_chart_data['vOTU'].iloc[other_idx_votu] + max(votu_chart_data['vOTU'][:-1]) * 0.02], 
        'k-', lw=1.5)

# Plot MAG chart with stacked bars
x_pos_mag = range(len(mag_chart_data))

# Plot total MAG (white bars)
ax2.bar(x_pos_mag, mag_chart_data['Host genome'], color='white', edgecolor='black', linewidth=0.5, label='Total MAG')

# Plot active MAG on top (black bars)
ax2.bar(x_pos_mag, mag_chart_data['MAG_is_active'], color='black', label='Active MAG')

# ax2.set_title('MAG Distribution (Top 20 + Other)', fontsize=14, fontweight='bold')
ax2.set_xlabel('Host Taxonomy Rank', fontsize=12)
ax2.set_ylabel('MAG Count', fontsize=12)
ax2.set_xticks(x_pos_mag)
ax2.set_xticklabels(mag_chart_data['highest_host_tax_rank'], rotation=45, ha='right')
ax2.legend()

# Add value labels on MAG bars (using actual values, not scaled)
for i, (scaled_val, actual_val) in enumerate(zip(mag_chart_data['Host genome'], mag_actual_values['Host genome'])):
    if mag_chart_data['highest_host_tax_rank'].iloc[i] == 'Other':
        # Position label above the scaled bar for "Other"
        ax2.text(i, scaled_val + max(mag_chart_data['Host genome'][:-1]) * 0.05, str(int(actual_val)), 
                ha='center', va='bottom', fontsize=9, fontweight='bold', color='red')
    else:
        ax2.text(i, scaled_val + 0.1, str(int(actual_val)), ha='center', va='bottom', fontsize=9)

# Add value labels for active MAG (using actual values)
for i, (scaled_total, scaled_active, actual_active) in enumerate(zip(
    mag_chart_data['Host genome'], mag_chart_data['MAG_is_active'], mag_actual_values['MAG_is_active'])):
    
    if actual_active > 0:
        if mag_chart_data['highest_host_tax_rank'].iloc[i] == 'Other':
            ax2.text(i, scaled_active/2, str(int(actual_active)), ha='center', va='center', 
                    color='pink', fontsize=9, fontweight='bold')
        else:
            ax2.text(i, scaled_active/2, str(int(actual_active)), ha='center', va='center', 
                    color='pink', fontsize=9, fontweight='bold')

# Add break indicator for "Other" in MAG chart
other_idx_mag = len(mag_chart_data) - 1
ax2.plot([other_idx_mag-0.4, other_idx_mag+0.4], 
        [mag_chart_data['Host genome'].iloc[other_idx_mag] + max(mag_chart_data['Host genome'][:-1]) * 0.02, 
         mag_chart_data['Host genome'].iloc[other_idx_mag] + max(mag_chart_data['Host genome'][:-1]) * 0.02], 
        'k-', lw=1.5)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_iphop_metadata.png', dpi=300)
plt.show()


# # Filter for active MAGs

# In[57]:


active_vOTUs_active_MAGs = (iphop_df
     .loc[
            (iphop_df['vOTU_is_active'] == True) 
          & (iphop_df['MAG_is_active'] == True)
         ]
)


# In[58]:


active_vOTUs_active_MAGs


# In[59]:


active_vOTUs_active_MAGs.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/active_vOTUs_active_MAGs_for_cytoscape.tsv', sep='\t', index=False)


# In[60]:


active_vOTUs_active_MAGs['highest_host_tax_rank'].nunique()


# In[61]:


active_vOTUs_active_MAGs.loc[active_vOTUs_active_MAGs['K'] == 'd__Archaea'].highest_host_tax_rank.nunique()


# In[62]:


active_vOTUs_active_MAGs.loc[active_vOTUs_active_MAGs['K'] == 'd__Bacteria'].highest_host_tax_rank.nunique()


# In[63]:


active_vOTUs_active_MAGs['P'].nunique()


# In[64]:


active_vOTUs_active_MAGs['S'].nunique()


# In[65]:


active_vOTUs_active_MAGs['G'].nunique()


# In[66]:


active_vOTUs_active_MAGs['Main method'].value_counts()


# In[70]:


active_vOTUs_active_MAGs


# In[71]:


active_vOTUs_active_MAGs['Main method'].value_counts().to_frame().reset_index()


# In[72]:


# How many predictions to the genus level?
len(active_vOTUs_active_MAGs.loc[active_vOTUs_active_MAGs['highest_host_tax_rank'].str.startswith('g__')]) / len(active_vOTUs_active_MAGs)


# # How many vOTUs have a host prediction out of how many are active?

# In[73]:


num_vOTU_hosts = (active_vOTUs_active_MAGs
.drop_duplicates(subset=['vOTU', 'Host genome'])
.groupby('highest_host_tax_rank')
.agg({
    'vOTU': 'nunique',
    'Host genome': 'nunique'
}
    )
).sort_values(by='vOTU', ascending=False)

num_vOTU_hosts


# In[74]:


import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

# Create the scatter plot
sns.scatterplot(data=num_vOTU_hosts, x='Host genome', y='vOTU', alpha=0.4)

# Calculate the regression line and correlation coefficient
slope, intercept, r_value, p_value, std_err = linregress(
    num_vOTU_hosts['Host genome'],  # Log-transform x-axis
    num_vOTU_hosts['vOTU']         # Log-transform y-axis
)

# Set log scale for both axes
plt.xscale('log')
plt.yscale('log')

r_squared = r_value**2

# Add R² text to the plot
plt.text(0.25, 0.95, f'R² = {r_squared:.3f}', transform=ax.transAxes, 
         fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.title('num vOTUs predicted to infect highest_host_tax_rank vs number of MAGs in highest_host_tax_rank')

# Display the plot
plt.show()


# # Methanogen viruses?

# In[62]:


active_vOTUs_active_MAGs.loc[active_vOTUs_active_MAGs['Host taxonomy'].str.contains('Methan')]


# In[63]:


iphop_df.loc[iphop_df['Host taxonomy'].str.contains('Methan')].vOTU.nunique()


# In[64]:


iphop_df.loc[(iphop_df['Host taxonomy'].str.contains('Methan')) & (iphop_df['vOTU_is_active'] == True)].vOTU.nunique()


# In[65]:


iphop_df.loc[(iphop_df['Host taxonomy'].str.contains('Methan')) & (iphop_df['vOTU_is_active'] == True)].vOTU.unique()


# In[66]:


iphop_df.loc[(iphop_df['Host taxonomy'].str.contains('Methan')) & (iphop_df['vOTU_is_active'] == True)].G.unique()


# # Check abundance correlation between host genomes and phage genomes, and genome length of phage genomes

# # Do more abundant hosts have more phages?

# In[67]:


host_taxonomy.columns = ['MAG', 'Host taxonomy']
active_MAGs['MAG'] = 'Add_' + active_MAGs['MAG']
active_MAGs


# In[68]:


active_MAGs = active_MAGs.merge(host_taxonomy, on='MAG', how='left')


# In[69]:


active_MAGs[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = active_MAGs['Host taxonomy'].str.split(';', expand=True)
active_MAGs['highest_host_tax_rank'] = active_MAGs.apply(get_highest_resolution, axis=1)


# In[70]:


active_MAGs = active_MAGs.drop(columns=['gene', 'K', 'P', 'C', 'O', 'F', 'G', 'S'])


# In[71]:


active_MAGs_melted = active_MAGs.melt(id_vars=['MAG', 'Host taxonomy', 'highest_host_tax_rank'], var_name='Sample', value_name='abundance')


# In[72]:


# get max and average
host_metaT_max = active_MAGs_melted.groupby('highest_host_tax_rank').agg({'abundance': 'max'}).reset_index().sort_values(by='abundance', ascending=False)
host_metaT_average = active_MAGs_melted.groupby('highest_host_tax_rank').agg({'abundance': 'mean'}).reset_index().sort_values(by='abundance', ascending=False)


# In[73]:


host_metaT_max = host_metaT_max.merge(num_vOTU_hosts, on='highest_host_tax_rank', how='inner')
host_metaT_average = host_metaT_average.merge(num_vOTU_hosts, on='highest_host_tax_rank', how='inner')


# In[74]:


import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

# Create the plot with regression line
ax = sns.regplot(data=host_metaT_max, x='abundance', y='vOTU', 
                 scatter=True, ci=None, line_kws={'color': 'red'})
# Calculate the regression line and correlation coefficient
slope, intercept, r_value, p_value, std_err = linregress(
    host_metaT_max['abundance'],
    host_metaT_max['vOTU']
)

# Set log scale for both axes
plt.xscale('log')
plt.yscale('log')

r_squared = r_value**2

# Add R² text to the plot
plt.text(0.05, 0.95, f'R² = {r_squared:.3f}', transform=ax.transAxes, 
         fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S4-vOTU_count_MAG_abundance_scatter.png', dpi=300)
# Display the plot
plt.show()


# In[75]:


import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

# Create the plot with regression line
ax = sns.regplot(data=host_metaT_max, x='abundance', y='Host genome', 
                 scatter=True, ci=None, line_kws={'color': 'red'})
# Calculate the regression line and correlation coefficient
slope, intercept, r_value, p_value, std_err = linregress(
    host_metaT_max['abundance'],
    host_metaT_max['Host genome']
)

# Set log scale for both axes
plt.xscale('log')
plt.yscale('log')

r_squared = r_value**2

# Add R² text to the plot
plt.text(0.05, 0.95, f'R² = {r_squared:.3f}', transform=ax.transAxes, 
         fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
# Display the plot
plt.show()


# In[76]:


import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

# Create the scatter plot
sns.scatterplot(data=host_metaT_average, x='abundance', y='vOTU', alpha=0.4)

# Calculate the regression line and correlation coefficient
slope, intercept, r_value, p_value, std_err = linregress(
    np.log10(host_metaT_average['abundance']),  # Log-transform x-axis
    np.log10(host_metaT_average['vOTU'])         # Log-transform y-axis
)

# Generate x-values for the regression line
x_values = np.logspace(
    np.log10(host_metaT_average['abundance'].min()), 
    np.log10(host_metaT_average['abundance'].max()), 
    100
)

# Calculate y-values for the regression line
y_values = 10 ** (slope * np.log10(x_values) + intercept)

# Plot the regression line
plt.plot(x_values, y_values, color='red', label=f'R² = {r_value**2:.2f}')

# Set log scale for both axes
plt.xscale('log')
plt.yscale('log')

# Add legend to display R²
plt.legend()

plt.title('num vOTUs predicted to infect highest_host_tax_rank vs average abundance')

# Display the plot
plt.show()


# In[ ]:




