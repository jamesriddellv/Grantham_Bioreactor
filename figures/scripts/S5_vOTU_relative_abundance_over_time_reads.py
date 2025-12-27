#!/usr/bin/env python
# coding: utf-8

# In[55]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


from print_versions import print_versions
print_versions(globals())


# In[56]:


df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
df.head()


# In[57]:


# Add gene key so we can filter the htseq table
gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/prodigal-gv/combined_manual_filtered_gene_ids.txt', sep='\t')
gff_df = gff_df[['vOTU', 'gene']]
gff_df['gene_position'] = gff_df['gene'].apply(lambda x: x.split('_')[1])
gff_df.head()


# In[58]:


# Read readcount data
data = pd.read_csv(
    '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED.tsv',
    sep='\t', header=None
)

# Assign column names
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

# Drop condensed tannin
data = data.drop(columns=['STM_0716_E_M_E029', 'STM_0716_E_M_E030', 'STM_0716_E_M_E031',
                       'STM_0716_E_M_E054', 'STM_0716_E_M_E055', 'STM_0716_E_M_E056',
                       'STM_0716_E_M_E066', 'STM_0716_E_M_E067', 'STM_0716_E_M_E068',
                       'STM_0716_E_M_E125', 'STM_0716_E_M_E126', 'STM_0716_E_M_E127'])

# Drop last two metadata rows
data = data.iloc[:-2, :]

# Add vOTU
data = data.merge(gff_df, on='gene', how='left')

# filter for only active vOTUs
data = data.loc[data['vOTU'].isin(df['vOTU'].unique())]



# In[59]:


vOTU_readcount = data.drop(columns=['gene', 'gene_position']).groupby('vOTU').sum()
vOTU_readcount


# In[60]:


df = vOTU_readcount


# In[61]:


df


# # Species accumulation curve

# In[62]:


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
plt.xticks(ticks=xtick_positions, labels=xtick_positions, fontsize=14)
plt.yticks(fontsize=14)

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

# In[63]:


# Split by catechin and unamended samples
replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[64]:


unamended_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'unamended') | (replicate_frame.day == 0)]['Sample'])
catechin_samples = list(replicate_frame.loc[(replicate_frame.treatment == 'catechin') | (replicate_frame.day == 0)]['Sample'])


# In[65]:


df_T_U = df_T.loc[df_T.index.isin(unamended_samples)]
df_T_C = df_T.loc[df_T.index.isin(catechin_samples)]


# In[66]:


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


# In[67]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/active_vOTUs_active_MAGs_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.sort_values(by=['vOTU', 'Confidence score']).drop_duplicates('vOTU')
# due to changes in taxonomy, we will change g__Lachnoclostridium to g__Clostridium
# iphop_df['G'] = iphop_df['G'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
# iphop_df['highest_host_tax_rank'] = iphop_df['highest_host_tax_rank'].str.replace('g__Lachnoclostridium', 'g__Clostridium')
iphop_df.head()


# In[68]:


# here we correctly handle situations where vOTUs might have had an iphop assignment, but wasn't an active MAG, therefore is called "no_assignment"
df = df.merge(iphop_df, on='vOTU', how='left')


# In[69]:


df.columns


# In[70]:


keep_cols = ['highest_host_tax_rank_appended', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']


# In[71]:


df['highest_host_tax_rank'] = df['highest_host_tax_rank'].fillna('no_assignment')
sum(df['highest_host_tax_rank'] != 'no_assignment') / len(df)


# In[72]:


abundant_vOTUs_list = ['contig_1022629',
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


# In[73]:


abundant_vOTUs = df.loc[df['vOTU'].isin(abundant_vOTUs_list)]


# In[74]:


def get_append_tax_rank(row):
    if row['vOTU'] in list(abundant_vOTUs['vOTU']):
        return row['highest_host_tax_rank'] + '__' + row['vOTU']
    else:
        return row['highest_host_tax_rank']


# In[75]:


df['highest_host_tax_rank_appended'] = df.apply(get_append_tax_rank, axis=1)
df


# In[76]:


# Sum counts in same rank
by_tax_rank = df.drop(columns=['highest_host_tax_rank']).groupby('highest_host_tax_rank_appended').sum().reset_index()[keep_cols]


# In[77]:


by_tax_rank_melted = by_tax_rank.melt(id_vars=['highest_host_tax_rank_appended'], var_name='Sample', value_name='sum_reads')
by_tax_rank_melted


# In[78]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')


# In[79]:


by_tax_rank_melted = by_tax_rank_melted.merge(replicate_frame, on='Sample', how='left')


# In[80]:


# add reads per sample
reads_per_sample = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv', usecols=['Sample', 'total reads (R1+R2)'])
by_tax_rank_melted = by_tax_rank_melted.merge(reads_per_sample, on='Sample', how='left')


# In[81]:


by_tax_rank_melted


# In[82]:


# Add prop abundance
by_sample_sum_abundance = by_tax_rank_melted.groupby(['Sample']).agg({'sum_reads': 'sum'}).reset_index().rename(columns={'sum_reads': 'sample_sum_reads'})
by_sample_sum_abundance


# In[83]:


by_tax_rank_melted = by_tax_rank_melted.merge(by_sample_sum_abundance, on='Sample', how='left')


# In[84]:


by_tax_rank_melted['prop_abundance'] = by_tax_rank_melted['sum_reads'] / by_tax_rank_melted['sample_sum_reads']
by_tax_rank_melted['prop_total'] = by_tax_rank_melted['sum_reads'] / by_tax_rank_melted['total reads (R1+R2)']


# In[85]:


by_tax_rank_melted


# In[86]:


by_tax_rank_melted['day'] = by_tax_rank_melted['day'].astype(int)


# In[87]:


# Step 1: Make sure 'taxon' is the correct column name

# Step 2: Group by treatment_day_replicate and taxon, calculate mean prop_abundance
# grouped = (
#     by_tax_rank_melted.groupby(['treatment_day_replicate', 'highest_host_tax_rank_appended'])['prop_abundance']
#     .mean()
#     .reset_index()
# )

# Step 3: Find max prop_abundance per taxon across all groups
max_abundance = by_tax_rank_melted.groupby('highest_host_tax_rank_appended')['prop_abundance'].max()

# Step 4: Identify taxa to keep (those with at least one group ≥ 5%) FROM GeTMM
taxa_to_keep = ['g__Clostridium', 'g__Clostridium__contig_507316',
       'g__Clostridium__contig_588865', 'g__Clostridium__contig_755337',
       'g__Fen-1039', 'g__JAGFXR01__contig_591846', 'g__PALSA-129',
       'g__Paludibacter__contig_724521', 'g__Pseudomonas_E',
       'g__Terracidiphilus', 'g__UBA1794', 'no_assignment',
       'no_assignment__contig_1022629', 'no_assignment__contig_1093732',
       'no_assignment__contig_709158', 'no_assignment__contig_713418',
       'no_assignment__contig_773239', 'no_assignment__contig_776767',
       'no_assignment__contig_861080']

print(taxa_to_keep)


# In[88]:


by_tax_rank_melted.columns


# In[89]:


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


# In[90]:


by_tax_rank_melted.columns


# In[91]:


by_tax_rank_melted_grouped = by_tax_rank_melted.groupby(['treatment', 'day', 'taxon_grouped']).agg(taxon_grouped_abundance=('sample_sum_reads', 'sum'), taxon_grouped_prop_total=('prop_total', 'mean')).reset_index()
by_tax_rank_melted_grouped


# In[92]:


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


# In[93]:


by_tax_rank_melted_grouped


# In[94]:


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


# In[95]:


# Pivot the dataframe: rows = day, columns = taxon, values = mean prop abundance
df_grouped = (
    by_tax_rank_melted_grouped
    .groupby(['treatment', 'day', 'taxon_grouped'])['taxon_grouped_prop_total']
    .sum()
    .reset_index()
).sort_values(by='taxon_grouped_prop_total', ascending=False)

# Choose a treatment to plot (e.g., 'unamended' or 'catechin')
treatment_to_plot = 'catechin'
df_treatment = df_grouped[df_grouped['treatment'] == treatment_to_plot]
# 1. Convert color dict to dataframe
color_df = pd.DataFrame(list(mag_color_dict.items()), columns=['taxon_grouped', 'color'])

# 2. Merge colors into df_treatment (long format)
df_treatment = df_treatment.merge(color_df, on='taxon_grouped', how='left')

# 3. Pivot to wide format (abundance table)
df_pivot = df_treatment.pivot(index='day', columns='taxon_grouped', values='taxon_grouped_prop_total')
df_pivot = df_pivot.fillna(0).sort_index()
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
plt.ylabel("proportion of active vOTU of total RNA", fontsize=14)
plt.xticks([0, 7, 14, 21, 35], fontsize=12)
# plt.yticks([0.0, 0.01, 0.02, 0.03], fontsize=12)
# plt.axvline(7, color='black')
# plt.title(f"{treatment_to_plot.capitalize()} Treatment", fontsize=16)
# plt.ylim(0, 1)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='Taxon')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(
    f"/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/stacked_area_{treatment_to_plot}_vOTU_metaT_TOTAL_simplified.pdf",
    dpi=300
)
plt.show()


# In[96]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__JAGFXR01'] * 100


# In[97]:


df.loc[df['highest_host_tax_rank'] == 'g__JAGFXR01']


# In[98]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__Pseudomonas_E'] * 100


# In[99]:


df_treatment.loc[df_treatment['taxon_grouped'] == 'g__Clostridium'] * 100


# In[100]:


df_pivot.sum(axis=1)


# # Plot individual phages from these >10%

# In[101]:


import matplotlib.colors as mcolors
treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

# Function to create a lighter color
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by mixing it with white.
    `amount` controls the mix ratio (0 = no change, 1 = white).
    """
    rgb = np.array(mcolors.to_rgb(color))  # Convert to RGB
    white = np.array([1, 1, 1])  # White color
    return mcolors.to_hex((1 - amount) * rgb + amount * white)  # Blend with white


def plot_vOTU_abundance_by_treatment(highest_host_tax_rank, df, with_individual=False, has_recombinase=False, is_active=False, with_replicate_points=True):
    
    if highest_host_tax_rank != []:
        df_subset = df.loc[
            (df['highest_host_tax_rank'].isin(highest_host_tax_rank)) & 
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()

    else:
        df_subset = df.loc[
            (df['treatment'].isin(['unamended', 'catechin']))
        ].copy()   

    if has_recombinase:
        df_subset = df_subset.loc[df['is_provirus_aggregated'] == True]

    if is_active:
        df_subset = df_subset.loc[df['is_active'] == True]
    
    sns.set_style("white")
    treatment_dict = {'unamended': '#7ACBC3', 'catechin': '#FC9B2D'}

    plt.figure(figsize=(8, 6))

    if with_individual:
        for (votu, treatment), group in df_subset.groupby(['vOTU', 'treatment']):
            # light_color = lighten_color(treatment_dict[treatment], amount=0.5)
            sns.lineplot(
                data=group, x='day', y='log2_abundance',
                color=treatment_dict[treatment], linewidth=0.8, alpha=0.5, errorbar=None, legend=False
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


# In[102]:


df


# In[103]:


df.columns


# In[104]:


df_melted = df[['vOTU', 'highest_host_tax_rank', 'STM_0716_E_M_E002', 'STM_0716_E_M_E003', 'STM_0716_E_M_E004',
       'STM_0716_E_M_E025', 'STM_0716_E_M_E027', 'STM_0716_E_M_E033',
       'STM_0716_E_M_E034', 'STM_0716_E_M_E035', 'STM_0716_E_M_E050',
       'STM_0716_E_M_E051', 'STM_0716_E_M_E052', 'STM_0716_E_M_E058',
       'STM_0716_E_M_E059', 'STM_0716_E_M_E060', 'STM_0716_E_M_E062',
       'STM_0716_E_M_E063', 'STM_0716_E_M_E064', 'STM_0716_E_M_E070',
       'STM_0716_E_M_E071', 'STM_0716_E_M_E072', 'STM_0716_E_M_E121',
       'STM_0716_E_M_E122', 'STM_0716_E_M_E123', 'STM_0716_E_M_E129',
       'STM_0716_E_M_E130', 'STM_0716_E_M_E131']].melt(id_vars=['vOTU', 'highest_host_tax_rank'], var_name='Sample', value_name='abundance')

df_melted = df_melted.merge(replicate_frame, on='Sample', how='left')
df_melted['log2_abundance'] = df_melted['abundance'].apply(lambda x: np.log2(x + 1))
df_melted.head()


# In[105]:


plot_vOTU_abundance_by_treatment(['g__JAGFXR01'], df_melted, with_individual=True)


# In[ ]:


plot_vOTU_abundance_by_treatment(['g__Paludibacter'], df_melted, with_individual=True)


# # label by catechin-enriched or catechin-sensitive

# In[ ]:


get_ipython().system('jupyter nbconvert --to script 004-vOTU-relative-abundance-over-time-reads.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S5_vOTU_relative_abundance_over_time_reads')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




