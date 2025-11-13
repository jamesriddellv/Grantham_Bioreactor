#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
df = df.set_index('vOTU')
df.head()


# # Species accumulation curve

# In[3]:


from itertools import permutations

# === Load your data ===
# === Transpose so columns are species ===
df_T = df.T  # Now rows = samples, columns = species

# === Updated species_accumulation function to return all data ===
def species_accumulation(df, permutations_count=100):
    species_acc = np.zeros((permutations_count, len(df)))

    for i in range(permutations_count):
        order = np.random.permutation(df.index)
        seen_species = set()
        counts = []

        for j, sample in enumerate(order):
            sample_species = df.loc[sample]
            new_species = sample_species[sample_species > 0].index
            seen_species.update(new_species)
            counts.append(len(seen_species))
        
        species_acc[i, :] = counts

    mean_acc = species_acc.mean(axis=0)
    min_acc = species_acc.min(axis=0)
    max_acc = species_acc.max(axis=0)
    
    return mean_acc, min_acc, max_acc

# === Compute curve data ===
mean_curve, min_curve, max_curve = species_accumulation(df_T)

# === Plot with min-max shaded range ===
sns.set(style="whitegrid")
x = np.arange(1, len(mean_curve) + 1)

ax = plt.figure(figsize=(6, 6))
sns.lineplot(x=x, y=mean_curve, label="Mean species richness", color='black')
plt.fill_between(x, min_curve, max_curve, alpha=0.3, label="Min–Max range", color='grey')

# Axis labels with larger font
plt.xlabel("Number of samples", fontsize=14)
plt.ylabel("vOTUs detected (cumulative)", fontsize=14)
xtick_positions = np.arange(1, 27, 5)
ytick_positions = np.arange(0, 901, 100)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(ticks=ytick_positions, labels=ytick_positions, fontsize=14)

# Start x-axis at 1
plt.xlim(left=1)

# Move legend to the right
# plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0, fontsize=12)
plt.legend().remove()
plt.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-B_specaccum_active.pdf", dpi=300)
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-B_specaccum_active.png", dpi=300)
plt.show()




# # Make for Catechin and Unamended

# In[4]:


# Split by catechin and unamended samples
replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[5]:


unamended_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'unamended') | (replicate_frame.day == 0)]['Sample'])
catechin_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'catechin') | (replicate_frame.day == 0)]['Sample'])


# In[6]:


df_T_U = df_T.loc[df_T.index.isin(unamended_samples)]
df_T_C = df_T.loc[df_T.index.isin(catechin_samples)]


# In[7]:


from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === Updated species_accumulation function to return all data ===
def species_accumulation(df, permutations_count=100):
    species_acc = np.zeros((permutations_count, len(df)))
    for i in range(permutations_count):
        order = np.random.permutation(df.index)
        seen_species = set()
        counts = []
        for j, sample in enumerate(order):
            sample_species = df.loc[sample]
            new_species = sample_species[sample_species > 0].index
            seen_species.update(new_species)
            counts.append(len(seen_species))
        
        species_acc[i, :] = counts
    mean_acc = species_acc.mean(axis=0)
    min_acc = species_acc.min(axis=0)
    max_acc = species_acc.max(axis=0)
    
    return mean_acc, min_acc, max_acc

# === Compute curve data for both treatments ===
mean_curve_U, min_curve_U, max_curve_U = species_accumulation(df_T_U)
mean_curve_C, min_curve_C, max_curve_C = species_accumulation(df_T_C)

# === Plot with min-max shaded range for both treatments ===
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))

# Unamended treatment
x_U = np.arange(1, len(mean_curve_U) + 1)
sns.lineplot(x=x_U, y=mean_curve_U, label="Unamended", color='#7bccc4', linewidth=2)
plt.fill_between(x_U, min_curve_U, max_curve_U, alpha=0.3, color='#7bccc4')

# Catechin treatment
x_C = np.arange(1, len(mean_curve_C) + 1)
sns.lineplot(x=x_C, y=mean_curve_C, label="Catechin", color='#fe9929', linewidth=2)
plt.fill_between(x_C, min_curve_C, max_curve_C, alpha=0.3, color='#fe9929')

# Axis labels with larger font
plt.xlabel("Number of samples", fontsize=14)
plt.ylabel("vOTUs detected (cumulative)", fontsize=14)

# Set x-ticks based on the maximum number of samples
max_samples = max(len(df_T_U), len(df_T_C))
xtick_positions = np.arange(1, max_samples + 1, 5)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(fontsize=14)

# Start x-axis at 1
plt.xlim(left=1)

# Add legend
plt.legend(fontsize=12, loc='lower right')
plt.grid(alpha=0.2)
plt.tight_layout()

# Uncomment to save figures
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1X_specaccum_active_by_treatment.pdf", dpi=300)
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1X_specaccum_active_by_treatment.png", dpi=300)
plt.show()


# In[8]:


from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))

# === Plot combined treatment curves ===
# Unamended treatment
x_U = np.arange(1, len(mean_curve_U) + 1)
sns.lineplot(x=x_U, y=mean_curve_U, label="Unamended", color='#7bccc4', linewidth=2)
plt.fill_between(x_U, min_curve_U, max_curve_U, alpha=0.3, color='#7bccc4')

# Catechin treatment
x_C = np.arange(1, len(mean_curve_C) + 1)
sns.lineplot(x=x_C, y=mean_curve_C, label="Catechin", color='#fe9929', linewidth=2)
plt.fill_between(x_C, min_curve_C, max_curve_C, alpha=0.3, color='#fe9929')

# === Plot overall curve ===
x_overall = np.arange(1, len(mean_curve) + 1)
sns.lineplot(x=x_overall, y=mean_curve, label="All samples", color='black', linewidth=2)
plt.fill_between(x_overall, min_curve, max_curve, alpha=0.2, color='grey')

# Axis labels with larger font
plt.xlabel("Number of samples", fontsize=14)
plt.ylabel("vOTUs detected (cumulative)", fontsize=14)

# Set x-ticks based on the maximum number of samples
max_samples = max(len(df_T_U), len(df_T_C), len(mean_curve))
xtick_positions = np.arange(1, max_samples + 1, 5)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(fontsize=14)

# Start x-axis at 1
xtick_positions = np.arange(1, 27, 5)
ytick_positions = np.arange(0, 901, 100)
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(ticks=ytick_positions, labels=ytick_positions, fontsize=14)

# Add legend
plt.legend(fontsize=12, loc='lower right')
plt.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1C_specaccum_active_total.pdf", dpi=300)
plt.show()


# In[9]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/active_vOTUs_active_MAGs_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.sort_values(by=['vOTU', 'Confidence score']).drop_duplicates('vOTU')
# due to changes in taxonomy, we will change g__Lachnoclostridium to g__Clostridium
# iphop_df['G'] = iphop_df['G'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
# iphop_df['highest_host_tax_rank'] = iphop_df['highest_host_tax_rank'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
iphop_df.head()


# In[10]:


# here we correctly handle situations where vOTUs might have had an iphop assignment, but wasn't an active MAG, therefore is called "no_assignment"
df = df.merge(iphop_df, on='vOTU', how='left')


# In[11]:


df


# In[12]:


df['highest_host_tax_rank'] = df['highest_host_tax_rank'].fillna('no_assignment')
sum(df['highest_host_tax_rank'] != 'no_assignment') / len(df)


# In[13]:


keep_cols = ['highest_host_tax_rank_appended', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']


# In[14]:


df.columns


# ### What are the most abundant vOTUs in each treatment?

# In[15]:


list(replicate_frame.loc[replicate_frame['treatment'] == 'unamended']['Sample'])


# In[16]:


# get sum_getmm but per treatment
catechin_samples = list(replicate_frame.loc[replicate_frame['treatment'] == 'catechin']['Sample'])
unamended_samples = list(replicate_frame.loc[replicate_frame['treatment'] == 'unamended']['Sample'])

df['max_catechin'] = df[catechin_samples].max(axis=1)
df['max_unamended'] = df[unamended_samples].max(axis=1)
df


# In[17]:


df.sort_values(by='max_catechin', ascending=False)


# In[18]:


df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/vOTU_metadata_getmm.tsv', sep='\t', index=False)


# ### Plot a rank abundance curve to find the cutoffs

# In[21]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Assuming you have your dataframe 'df' with the columns
# df = pd.read_csv('your_data.csv')  # Load your data here

cutoff_dict = {
    'Max Catechin - Rank Abundance': 25000,
    'Max Unamended - Rank Abundance': 1000
}

def find_elbow_point(values):
    """
    Find the elbow point in a rank abundance curve using the "elbow method"
    This finds the point of maximum curvature
    """
    ranks = np.arange(1, len(values) + 1)
    
    # Calculate the distances from each point to the line connecting first and last points
    # This is a common method for finding the elbow
    first_point = np.array([ranks[0], values.iloc[0]])
    last_point = np.array([ranks[-1], values.iloc[-1]])
    
    distances = []
    for i in range(len(ranks)):
        point = np.array([ranks[i], values.iloc[i]])
        # Calculate perpendicular distance from point to line
        distance = np.abs(np.cross(last_point - first_point, first_point - point)) / np.linalg.norm(last_point - first_point)
        distances.append(distance)
    
    # The elbow is at the point with maximum distance
    elbow_index = np.argmax(distances)
    return ranks[elbow_index], values.iloc[elbow_index]

def plot_rank_abundance(values, title, ax):
    """
    Create a rank abundance curve with elbow point marked
    """
    # Remove NaN values and sort in descending order
    clean_values = values.dropna().sort_values(ascending=False)
    ranks = np.arange(1, len(clean_values) + 1)
    
    # Plot the actual data
    ax.plot(ranks, clean_values, 'o-', linewidth=2, markersize=4, label='Rank curve', color='blue')
    
    # Find and mark the elbow point
    elbow_rank, elbow_value = find_elbow_point(clean_values)
    ax.axvline(x=elbow_rank, color='red', linestyle='--', linewidth=2, alpha=0.8, 
               label=f'Elbow at rank {elbow_rank:.0f}')
    ax.plot(elbow_rank, elbow_value, 'ro', markersize=8, label=f'Elbow point')
    
    ax.set_xlabel('Rank')
    ax.set_ylabel('Value')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Add some statistics to the plot
    mean_val = clean_values.mean()
    median_val = clean_values.median()
    cutoff = cutoff_dict[title]
    ax.axhline(y=mean_val, color='green', linestyle=':', alpha=0.7, label=f'Mean: {mean_val:.2f}')
    ax.axhline(y=median_val, color='orange', linestyle=':', alpha=0.7, label=f'Median: {median_val:.2f}')
    ax.axhline(y=cutoff, color='red', linestyle='-', alpha=0.9, label='cutoff for individual plotting')
    
    ax.legend()
    
    return clean_values, elbow_rank

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot rank abundance curves
catechin_values, catechin_elbow = plot_rank_abundance(df['max_catechin'], 'Max Catechin - Rank Abundance', ax1)
unamended_values, unamended_elbow = plot_rank_abundance(df['max_unamended'], 'Max Unamended - Rank Abundance', ax2)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/SX_rank_abundance_cutoff_individual_plotting.pdf', dpi=300)
plt.show()


# Based on observations here, let's individually plot vOTUs where max catechin is >25000 and max unamended is >1000

# In[19]:


abundant_vOTUs = df.loc[(df['max_catechin'] > 25000) | (df['max_unamended'] > 1000)]


# In[23]:


list(abundant_vOTUs['vOTU'])


# In[33]:


abundant_vOTUs_host = list(abundant_vOTUs['highest_host_tax_rank'] + '__' + abundant_vOTUs['vOTU'])


# In[34]:


abundant_vOTUs_host


# In[35]:


def get_append_tax_rank(row):
    if row['vOTU'] in list(abundant_vOTUs['vOTU']):
        return row['highest_host_tax_rank'] + '__' + row['vOTU']
    else:
        return row['highest_host_tax_rank']


# In[36]:


df['highest_host_tax_rank_appended'] = df.apply(get_append_tax_rank, axis=1)
df


# In[37]:


# Sum getmm values of vOTUs in the same highest host tax rank
by_tax_rank = df.drop(columns=['highest_host_tax_rank']).groupby('highest_host_tax_rank_appended').sum().reset_index()[keep_cols]


# In[38]:


by_tax_rank_melted = by_tax_rank.melt(id_vars=['highest_host_tax_rank_appended'], var_name='Sample', value_name='abundance')
by_tax_rank_melted


# In[39]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# In[40]:


by_tax_rank_melted = by_tax_rank_melted.merge(replicate_frame, on='Sample', how='left')


# In[41]:


by_tax_rank_melted


# In[42]:


# Take the mean across replicates, since there are fewer replicates of unamended day 14
by_tax_rank_melted = by_tax_rank_melted.groupby(['highest_host_tax_rank_appended', 'treatment', 'day']).agg({'abundance': 'mean'}).reset_index()


# In[43]:


by_tax_rank_melted


# In[44]:


# Add prop abundance
by_sample_sum_abundance = by_tax_rank_melted.groupby(['treatment', 'day']).agg({'abundance': 'sum'}).reset_index().rename(columns={'abundance': 'total_getmm'})
by_sample_sum_abundance


# In[45]:


by_tax_rank_melted = by_tax_rank_melted.merge(by_sample_sum_abundance, on=['treatment', 'day'], how='left')


# In[46]:


by_tax_rank_melted['prop_abundance'] = by_tax_rank_melted['abundance'] / by_tax_rank_melted['total_getmm']


# In[47]:


by_tax_rank_melted


# In[48]:


by_tax_rank_melted['day'] = by_tax_rank_melted['day'].astype(int)


# In[49]:


# Step 1: Make sure 'taxon' is the correct column name

# Step 2: Group by treatment_day_replicate and taxon, calculate mean prop_abundance
# grouped = (
#     by_tax_rank_melted.groupby(['treatment_day_replicate', 'highest_host_tax_rank_appended'])['prop_abundance']
#     .mean()
#     .reset_index()
# )

# Step 3: Find max prop_abundance per taxon across all groups
max_abundance = by_tax_rank_melted.groupby('highest_host_tax_rank_appended')['prop_abundance'].max()

# Step 4: Identify taxa to keep (those with at least one group ≥ 5%)
taxa_to_keep = max_abundance[max_abundance >= 0.04].index

print(taxa_to_keep)


# In[50]:


set(taxa_to_keep) - set(abundant_vOTUs_host)


# In[51]:


taxa_to_keep = ['Other',
 'g__Clostridium',
 'g__Clostridium__contig_507316',
 'g__Clostridium__contig_588865',
 'g__Clostridium__contig_755337',
 'g__Fen-1039',
 'g__Fen-1308',
 'g__JAGFXR01',
 'g__JAGFXR01__contig_591846',
 'g__PALSA-129',
 'g__Paludibacter',
 'g__Paludibacter__contig_724521',
 'g__Pseudomonas_E',
 'g__Terracidiphilus',
 'g__UBA1794',
 'no_assignment',
 'no_assignment__contig_1022629',
 'no_assignment__contig_1093732',
 'no_assignment__contig_709158',
 'no_assignment__contig_713418',
 'no_assignment__contig_773239',
 'no_assignment__contig_776767',
 'no_assignment__contig_861080']

# Step 5: Create a new column in original df with taxon or "other"
by_tax_rank_melted['taxon_grouped'] = by_tax_rank_melted['highest_host_tax_rank_appended'].where(by_tax_rank_melted['highest_host_tax_rank_appended'].isin(taxa_to_keep), 'Other')


# In[52]:


by_tax_rank_melted.columns


# In[53]:


by_tax_rank_melted


# In[54]:


by_tax_rank_melted_grouped = by_tax_rank_melted.groupby(['treatment', 'day', 'taxon_grouped']).agg(taxon_grouped_abundance=('abundance', 'sum'), taxon_grouped_prop_abundance=('prop_abundance', 'sum')).reset_index()
by_tax_rank_melted_grouped


# In[55]:


# Step 1: Filter unamended day 0 rows
day0_unamended = by_tax_rank_melted_grouped[
    (by_tax_rank_melted_grouped['treatment'] == 'unamended') &
    (by_tax_rank_melted_grouped['day'] == 0)
].copy()

# Step 2: Modify treatment to 'catechin'
day0_unamended['treatment'] = 'catechin'

# Optional: also update identifiers like 'sample_treatment_day' or 'treatment_day_replicate' if needed
# For example:
# day0_unamended['sample_treatment_day'] = day0_unamended['sample_treatment_day'].str.replace('unamended', 'catechin')
# day0_unamended['treatment_day_replicate'] = day0_unamended['treatment_day_replicate'].str.replace('unamended', 'catechin')

# Step 3: Append to original dataframe
by_tax_rank_melted_grouped = pd.concat([by_tax_rank_melted_grouped, day0_unamended], ignore_index=True)


# In[56]:


by_tax_rank_melted_grouped


# In[57]:


mag_color_dict = {
    'Other': 'black',
    
    # g__Clostridium group (red shades) - maintaining existing colors
    'g__Clostridium': '#9E0B0F',                    # Original sangria red
    'g__Clostridium__contig_507316': '#DC143C',     # Crimson (from original)
    'g__Clostridium__contig_588865': '#FF4500',     # Orange red (new - similar shade)
    'g__Clostridium__contig_755337': '#B22222',     # Fire brick red (from original)
    
    # g__Fen-1039 group (orange shades)
    'g__Fen-1039': '#FF8C00',                       # Original dark orange
    
    # g__JAGFXR01 group (purple shades) - maintaining existing colors
    'g__JAGFXR01': '#845EC2',                       # Original purple
    'g__JAGFXR01__contig_591846': '#5C4188',        # Darker purple (from original)
    
    # g__PALSA-129 group (blue shades) - new from second list
    'g__PALSA-129': '#004E89',                      # Original deep blue
    
    # g__Methanosarcina - from first list
    'g__Methanosarcina': '#228B22',                 # Forest green (from original)
    
    # g__Paludibacter group (bright red shades) - maintaining existing colors
    'g__Paludibacter': '#FF1744',                   # Original bright red
    'g__Paludibacter__contig_724521': '#E91E63',    # Pink-red (from original)
    
    # g__Pseudomonas_E group (gold shades)
    'g__Pseudomonas_E': '#FFD700',                  # Original gold
    
    # g__Terracidiphilus group (green shades)
    'g__Terracidiphilus': '#32CD32',                # Original lime green
    
    # g__UBA1794 group (pink shades)
    'g__UBA1794': '#FF69B4',                        # Original hot pink
    
    # no_assignment group (gray shades) - maintaining existing colors
    'no_assignment__contig_1022629': '#C0C0C0',     # Silver (from original)
    'no_assignment__contig_1093732': '#D3D3D3',     # Light gray (from original)
    'no_assignment__contig_709158': '#696969',      # Dim gray (from original)
    'no_assignment__contig_713418': '#778899',      # Light slate gray (from original)
    'no_assignment__contig_773239': '#708090',      # Slate gray (from original)
    'no_assignment__contig_776767': '#556B2F',      # Dark olive green (from original)
    'no_assignment__contig_861080': '#2F4F4F',      # Dark slate gray (from original)
    'no_assignment': '#808080',                     # Medium gray
}


# In[58]:


list(mag_color_dict.keys())


# In[59]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

# Set yticks with custom formatting
import matplotlib.ticker as ticker

# Define the tick values and labels
yticks_values = [0, 5000, 10000, 15000, 20000]
yticks_labels = ['0', '5k', '10k', '15k', '20k']

# Set the yticks
plt.yticks(yticks_values, yticks_labels, fontsize=12)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[60]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

# Set yticks with custom formatting
import matplotlib.ticker as ticker

# Define the tick values and labels
yticks_values = [0, 5000, 10000, 15000, 20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000]
yticks_labels = ['0', '5k', '10k', '15k', '20k', '40k', '60k', '80k', '100k', '120k', '140k', '160k']

# Set the yticks
plt.yticks(yticks_values, yticks_labels, fontsize=12)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[61]:


## Create a separate figure just for the legend
fig_legend = plt.figure(figsize=(8, 6))
# Create proxy artists (patches) for the legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=mag_color_dict.get(col, 'gray'), 
                        label=col, alpha=0.7) 
                  for col in df_pivot.columns]
# Reverse the order
legend_elements = legend_elements[::-1]
# Create the legend
legend = fig_legend.legend(handles=legend_elements, 
                          loc='center',
                          ncol=1,  # Adjust number of columns as needed
                          fontsize=20)
                          # frameon=True)
                          # fancybox=True,
                          # shadow=True)
# Remove axes
fig_legend.gca().set_axis_off()
# Save the legend
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_legend_only_simplified.pdf",
    dpi=300,
    bbox_inches='tight'
)
plt.show()


# In[62]:


df.loc[df['highest_host_tax_rank'] == 'g__JAGFXR01']


# In[63]:


102652.241407 / 750.845462


# In[64]:


by_tax_rank_melted_grouped


# In[65]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'unamended'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[66]:


# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = by_tax_rank_melted_grouped[by_tax_rank_melted_grouped['treatment'] == treatment_to_plot]
# df_treatment = df_treatment.sort_values(by='taxon_grouped_prop_abundance', ascending=False)

# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_abundance')
# df_pivot = df_pivot.fillna(0).sort_index()
# df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)

# 4. Get colors in the exact column order of df_pivot
colors = [mag_color_dict.get(col, 'gray') for col in df_pivot.columns]

# 5. Plot
sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))
plt.stackplot(
    df_pivot.index,
    df_pivot.T.values,
    labels=df_pivot.columns,
    colors=colors,
    alpha=0.7,
    edgecolor='none'
)
plt.xlabel("Day", fontsize=14)
plt.ylabel("vOTU relative abundance by predicted host genus (geTMM)", fontsize=12)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)

plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-C_vOTU_stacked_area_prop_{treatment_to_plot}.pdf",
    dpi=300
)
plt.show()


# In[ ]:





# In[ ]:




