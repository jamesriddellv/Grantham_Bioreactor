#!/usr/bin/env python
# coding: utf-8

# In[1]:


# plot clostridium, UBA-1794, and figure out if there are different MAGs enriched at the different timepoints.
# See if these clusters align with the vOTUs being different, and see if the predictions line up as well.


# In[2]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[3]:


metaT = pd.read_excel('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/MAG_abundance_table.xlsx')
metaT.head()


# In[4]:


# Drop Condensed Tannin samples
metaT = metaT.drop(columns=['STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                        'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                        'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                        'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127'])


# In[5]:


# drop zeroes
metaT = metaT[metaT.iloc[:, 1:-1].sum(axis=1) > 0]

# Reshape to long format
long_df = metaT.melt(id_vars=["MAG", "GTDB"], 
                  var_name="Sample", 
                  value_name="abundance")
long_df


# In[6]:


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
    'STM_0716_E_M_E121_unamended_7': 1,
       'STM_0716_E_M_E122_unamended_7': 2,
    'STM_0716_E_M_E123_unamended_7': 3,
       'STM_0716_E_M_E125_CT_7': 1,
    'STM_0716_E_M_E126_CT_7': 2,
       'STM_0716_E_M_E127_CT_7': 3,
    'STM_0716_E_M_E129_catechin_7': 1,
       'STM_0716_E_M_E130_catechin_7': 2,
    'STM_0716_E_M_E131_catechin_7': 3
}.items())

replicate_frame['Sample'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[0])
replicate_frame['treatment'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[1])
replicate_frame['day'] = replicate_frame[0].apply(lambda x: x.rsplit('_', 2)[2])
replicate_frame.columns = ['sample_treatment_day', 'replicate', 'Sample', 'treatment', 'day']
replicate_frame['treatment_day_replicate'] = replicate_frame['treatment'] + '_' + replicate_frame['day'].astype(str) + '_' + replicate_frame['replicate'].astype(str)


# In[7]:


melted_df = long_df.merge(replicate_frame, on='Sample', how='left')
melted_df['replicate'] = melted_df['replicate'].astype(str)
melted_df.head()


# In[10]:


melted_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/getmms_REV_active_MAGs_melted.tsv', sep='\t', index=False)


# In[11]:


melted_df.MAG.nunique()


# In[12]:


active_MAGs = melted_df


# In[13]:


active_MAGs['day'] = active_MAGs['day'].astype(int)


# In[14]:


active_MAGs


# In[15]:


active_MAGs[['K', 'P', 'C', 'O', 'F', 'G', 'S']] = active_MAGs['GTDB'].str.split(';', expand=True)


# In[16]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

active_MAGs['highest_host_tax_rank'] = active_MAGs.apply(get_highest_resolution, axis=1)


# In[17]:


active_MAGs.loc[active_MAGs['highest_host_tax_rank'] == 'g__Methanosarcina'].MAG.unique()


# In[18]:


min_value = active_MAGs['abundance'].min()
shift_value = 1 - min_value
active_MAGs['log2_abundance'] = active_MAGs['abundance'].apply(lambda x: np.log2(x+shift_value))


# In[19]:


active_MAGs.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaT/active_MAGs_abundance_with_metadata.tsv', sep='\t', index=False)


# In[20]:


active_MAGs.head()


# In[21]:


data_mtx = active_MAGs.pivot_table(values='abundance', index='Sample', columns='MAG')


# In[24]:


# Read metadata
metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")

# Simplify treatment labels
metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('y')[1]))


# In[25]:


metadata = metadata.loc[metadata['treatment'] != 'T']


# In[26]:


# Plot MAG diversity statistics
from skbio.diversity.alpha import simpson
metadata['inverse_simpson'] = data_mtx.apply(lambda x: simpson(x), axis=1).values


# In[27]:


# Plot MAG diversity statistics
from skbio.diversity.alpha import shannon
metadata['shannon'] = data_mtx.apply(lambda x: shannon(x), axis=1).values


# In[28]:


metadata


# In[30]:


# Plot individual points and a mean line per treatment
plt.figure(figsize=(6, 6))

# Points (individual replicates)
sns.scatterplot(
    data=metadata,
    x='timepoint',
    y='shannon',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    alpha=1
)

# Mean line across timepoints per treatment
sns.lineplot(
    data=metadata,
    x='timepoint',
    y='shannon',
    hue='treatment',
    palette={"C": "#FC9B2D", "U": "#7ACBC3"},
    estimator='mean',
    errorbar='ci',
    linewidth=2,
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
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_MAG_alpha_diversity.png', dpi=300)
plt.show()



# In[ ]:


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
plt.ylabel("Simpson Diversity Index (1-D)", fontsize=14)
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
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/simpson_diversity_MAG.pdf', dpi=300)
plt.show()


# In[ ]:


# Function to create a lighter color
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by mixing it with white.
    `amount` controls the mix ratio (0 = no change, 1 = white).
    """
    rgb = np.array(mcolors.to_rgb(color))  # Convert to RGB
    white = np.array([1, 1, 1])  # White color
    return mcolors.to_hex((1 - amount) * rgb + amount * white)  # Blend with white


# In[ ]:


def plot_MAG_abundance_by_treatment(highest_host_tax_rank, with_individual=False, with_replicate_points=True, df=active_MAGs):
    
    if highest_host_tax_rank != []:
        df_subset = df.loc[
            (df['highest_host_tax_rank'].isin(highest_host_tax_rank)) & 
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()

    else:
        df_subset = df.loc[
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()   

    sns.set_style("white")
    treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

    plt.figure(figsize=(8, 6))

    if with_individual:
        for (votu, treatment), group in df_subset.groupby(['MAG', 'treatment']):
            light_color = lighten_color(treatment_dict[treatment], amount=0.5)
            sns.lineplot(
                data=group, x='day', y='log2_abundance',
                color=light_color, linewidth=0.8, alpha=0.5, errorbar=None, legend=False
            )


    # Plot mean lines per treatment
    for treatment, color in treatment_dict.items():
        df_treat = df_subset[df_subset['treatment'] == treatment]
        sns.lineplot(
            data=df_treat, x='day', y='log2_abundance',
            color=color, linewidth=2.5, errorbar=('ci', 95), label=treatment
        )


    if with_replicate_points and 'replicate' in df_subset.columns:
        # Plot replicate means as points
        mean_by_replicate = (
            df_subset.groupby(['treatment', 'day', 'replicate'])['log2_abundance']
            .mean()
            .reset_index()
        )
        for treatment, color in treatment_dict.items():
            subset = mean_by_replicate[mean_by_replicate['treatment'] == treatment]
            sns.scatterplot(
                data=subset, x='day', y='log2_abundance',
                color=color, s=60, edgecolor='black', linewidth=0.5,
                alpha=0.7, label=None
            )
    
    plt.xlabel("Day", fontsize=12)
    plt.ylabel("log2(Abundance)", fontsize=12)
    # plt.title(f'vOTU Abundance Over Time ({highest_host_tax_rank})', fontsize=14)
    plt.legend(title='Treatment')
    plt.legend().remove()
    plt.tight_layout()
    plt.show()


# In[ ]:


import seaborn as sns
import matplotlib.pyplot as plt

def plot_MAG_abundance_by_treatment(highest_host_tax_rank, with_individual=False, with_replicate_points=True, df=active_MAGs, shareyaxis=True):
    # Filter the dataframe based on host tax rank and treatment
    if highest_host_tax_rank:
        df_subset = df.loc[
            (df['highest_host_tax_rank'].isin(highest_host_tax_rank)) &
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()
    else:
        df_subset = df.loc[df['treatment'].isin(['unamended', 'catechin'])].copy()

    # Safety check: ensure we have more than one MAG
    if 'MAG' not in df_subset.columns or df_subset['MAG'].nunique() <= 1:
        print("Not enough unique MAGs to facet.")
        return

    sns.set_style("white")
    treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

    # Optional: precompute replicate means if needed
    if with_replicate_points and 'replicate' in df_subset.columns:
        replicate_means = (
            df_subset.groupby(['MAG', 'treatment', 'day', 'replicate'])['log2_abundance']
            .mean()
            .reset_index()
        )
    else:
        replicate_means = None

    # Create FacetGrid
    g = sns.FacetGrid(df_subset, col='MAG', col_wrap=3, sharey=shareyaxis, height=4, aspect=1.3)

    def facet_plot(data, **kwargs):
        mag = data['MAG'].iloc[0]
        for treatment in data['treatment'].unique():
            subset = data[data['treatment'] == treatment]
            # Mean line
            sns.lineplot(
                data=subset, x='day', y='log2_abundance',
                color=treatment_dict[treatment], linewidth=2.0, errorbar=('ci', 95), label=treatment
            )
            # Individual lines
            if with_individual:
                for replicate_id, grp in subset.groupby('replicate'):
                    light_color = lighten_color(treatment_dict[treatment], amount=0.5)
                    sns.lineplot(
                        data=grp, x='day', y='log2_abundance',
                        color=light_color, linewidth=0.8, alpha=0.5, errorbar=None, legend=False
                    )
            # Replicate points
            if replicate_means is not None:
                points = replicate_means[(replicate_means['MAG'] == mag) & (replicate_means['treatment'] == treatment)]
                sns.scatterplot(
                    data=points, x='day', y='log2_abundance',
                    color=treatment_dict[treatment], s=60, edgecolor='black',
                    linewidth=0.5, alpha=0.7, legend=False
                )

    g.map_dataframe(facet_plot)

    g.set_axis_labels("Day", "log2(Abundance)")
    g.set_titles(col_template="{col_name}")
    g.tight_layout()
    plt.show()


# In[ ]:


plot_MAG_abundance_by_treatment(['g__UBA1794'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Methanobacterium_B'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Clostridium'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__JAGFXR01'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Polaromonas'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Solibacter'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Methanosarcina'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__LD21'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__LD21'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Methanobacterium_A'])


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Methanobacterium_B'])


# In[ ]:


melted_df.loc[melted_df['MAG'] == 'STM_0716_E_M_E050_E054_E058_D_bin.80']


# In[ ]:


plot_MAG_abundance_by_treatment(['g__Terracidiphilus'])


# In[ ]:





# In[18]:


metaG = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/meancov/metaG/MAGs/metaG_MAG_mean_97_75_75.tsv', sep='\t')
metaG.head()


# In[19]:


# Drop zero rows
present_MAGs = metaG[metaG.iloc[:, 1:].sum(axis=1) > 0]

# Read in sample metadata
sample_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaG_sample_metadata.csv')

# Clean up column names to match sample_metadata
present_MAGs.columns = [col.replace("_97FILTERED_SORTED Mean", "").replace("_MG", "") for col in present_MAGs.columns]

# Melt the DataFrame to long format
present_MAGs = present_MAGs.melt(id_vars=['Genome'], var_name='Sample', value_name='abundance')

# Add metadata to the melted DataFrame
sample_metadata['sample_treatment_day'] = sample_metadata['Sample'] + "_" + sample_metadata['treatment'] + "_" + sample_metadata['timepoint'].str.replace("day", "")
present_MAGs = present_MAGs.merge(sample_metadata, left_on='Sample', right_on='Sample', how='left')
present_MAGs['day'] = present_MAGs['timepoint'].str.replace("day", "").astype(int)

present_MAGs = present_MAGs[['Genome', 'abundance', 'Sample', 'sample_treatment_day', 'treatment', 'day']]

# add replicates
sample_replicate = active_MAGs[['Sample', 'replicate']].drop_duplicates()
present_MAGs = present_MAGs.merge(sample_replicate, on='Sample', how='left').fillna('4')
present_MAGs = present_MAGs.rename(columns={'Genome': 'MAG'})
# Display the final DataFrame
present_MAGs


# In[20]:


present_MAGs = present_MAGs.merge(active_MAGs[['MAG', 'GTDB', 'K', 'P', 'C', 'O', 'F', 'G', 'S']].drop_duplicates(), on='MAG', how='left')
present_MAGs


# In[21]:


# Define the order of columns from highest to lowest resolution
taxonomy_columns = ['G', 'F', 'O', 'C', 'P']

# Create a new column with the highest resolution taxonomic value
def get_highest_resolution(row):
    for col in taxonomy_columns:
        if row[col] != f'{col.lower()}__':  # Check if the value is not the default prefix
            return row[col]
    return np.nan  # If no valid value is found, return NaN

present_MAGs['highest_host_tax_rank'] = present_MAGs.apply(get_highest_resolution, axis=1)


# In[22]:


present_MAGs.head()


# In[23]:


min_value = present_MAGs['abundance'].min()
shift_value = 1 - min_value
present_MAGs['log2_abundance'] = present_MAGs['abundance'].apply(lambda x: np.log2(x+shift_value))


# In[24]:


present_MAGs.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/present_MAGs_abundance_with_metadata.tsv', sep='\t', index=False)


# In[25]:


plot_MAG_abundance_by_treatment(['g__JAGFXR01'], df=present_MAGs)


# In[26]:


plot_MAG_abundance_by_treatment(['g__Clostridium'], df=present_MAGs)


# In[27]:


plot_MAG_abundance_by_treatment(['g__Pelorhabdus'], df=present_MAGs)


# In[28]:


plot_MAG_abundance_by_treatment(['g__Pseudomonas_E'], df=present_MAGs)


# In[32]:


import seaborn as sns
import matplotlib.pyplot as plt

def plot_average_MAG_abundance_by_treatment(highest_host_tax_rank, show_individual_MAGs=False, df=active_MAGs):
    # Filter the dataframe
    if highest_host_tax_rank:
        df_subset = df.loc[
            (df['highest_host_tax_rank'].isin(highest_host_tax_rank)) &
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()
    else:
        df_subset = df.loc[df['treatment'].isin(['unamended', 'catechin'])].copy()

    sns.set_style("white")
    treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

    plt.figure(figsize=(8, 6))

    # Optionally plot all individual MAGs as faint lines
    if show_individual_MAGs:
        for mag in df_subset['MAG'].unique():
            for treatment in ['unamended', 'catechin']:
                subset = df_subset[(df_subset['MAG'] == mag) & (df_subset['treatment'] == treatment)]
                if not subset.empty:
                    light_color = lighten_color(treatment_dict[treatment], amount=0.5)
                    sns.lineplot(
                        data=subset, x='day', y='log2_abundance',
                        color=light_color, linewidth=0.8, alpha=0.4, errorbar=None
                    )

    # Compute MAG-wise means per day & treatment
    mag_means = (
        df_subset.groupby(['treatment', 'day', 'MAG'])['log2_abundance']
        .mean()
        .reset_index()
    )

    # Now compute average over MAGs, retaining variation between MAGs
    for treatment in ['unamended', 'catechin']:
        subset = mag_means[mag_means['treatment'] == treatment]
        sns.lineplot(
            data=subset, x='day', y='log2_abundance',
            color=treatment_dict[treatment], linewidth=2.5, errorbar=('ci', 95), label=treatment
        )

    plt.xlabel("Day", fontsize=12)
    plt.ylabel("Mean log2(Abundance) across MAGs", fontsize=12)
    plt.legend(title='Treatment')
    plt.tight_layout()
    plt.show()


# In[33]:


plot_average_MAG_abundance_by_treatment(['g__Pseudomonas_E'], df=present_MAGs)


# In[35]:


plot_average_MAG_abundance_by_treatment(['g__Pseudomonas_E'], df=active_MAGs)


# In[84]:


sample_metadata


# In[85]:


present_MAGs.loc[present_MAGs['G'] == 'g__Pseudomonas_E'].sort_values(by='log2_abundance', ascending=False).drop_duplicates()


# In[ ]:





# # How many reads mapped to MAGs?

# In[48]:


count = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/meancov/metaG/MAGs/count/metaG_MAG_count_97_75.tsv", sep='\t')
count.head()


# In[49]:


count.columns = [j.split('MG')[0].strip('_') for j in [i.split('_97')[0] for i in count.columns]]
count.head()


# In[50]:


metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaG_sample_metadata.csv')
metadata


# In[51]:


reads_mapped = count.iloc[:,1:].sum(axis=0).reset_index()
reads_mapped.columns = ['Sample', 'reads mapped to MAGs']
reads_mapped


# In[53]:


reads_mapped = metadata.merge(reads_mapped, on='Sample')


# In[54]:


reads_mapped


# In[57]:


reads_mapped['prop_mapped_to_MAGs'] = reads_mapped['reads mapped to MAGs'] / reads_mapped['total reads (R1+R2)']


# In[58]:


reads_mapped


# In[ ]:




