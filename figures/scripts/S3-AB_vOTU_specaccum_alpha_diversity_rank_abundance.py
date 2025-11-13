#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import shannon
from skbio.diversity.alpha import observed_otus
import seaborn as sns


# # This script investigates what species accumulation and alpha diversity would look like under different activity cutoffs

# In[3]:


df_1_read_mapped = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_vOTUs_wide.tsv', sep='\t')
df_prop10 = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_prop10_vOTUs_wide.tsv', sep='\t')
df_1_gene_mapped_per_10kb = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1_gene_per_10kb_vOTUs_wide.tsv', sep='\t')
df_5_reads_mapped = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')


# # Species accumulation

# In[4]:


from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def species_accumulation(df, permutations_count=100):
    """
    Calculate species accumulation curves with permutation-based confidence intervals.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with samples as rows and species/vOTUs as columns
    permutations_count : int
        Number of permutations to run for confidence intervals
        
    Returns:
    --------
    tuple : (mean_curve, min_curve, max_curve)
        Mean, minimum, and maximum accumulation curves
    """
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

def prepare_treatment_data_for_single_df(df, treatment_cols_dict):
    """
    Prepare treatment data for a single DataFrame.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with coverage data
    treatment_cols_dict : dict
        Dictionary mapping treatment names to column lists
        
    Returns:
    --------
    dict : Dictionary with treatment DataFrames ready for accumulation analysis
    """
    treatment_dfs = {}
    
    for treatment_name, columns in treatment_cols_dict.items():
        # Filter columns that exist in the DataFrame
        common_cols = [col for col in columns if col in df.columns]
        
        # Create DataFrame for this treatment
        df_treatment = df[common_cols].T  # Transpose for species accumulation function
        treatment_dfs[treatment_name] = df_treatment
    
    return treatment_dfs

def plot_single_species_accumulation(ax, treatment_curves, title, colors=None):
    """
    Plot species accumulation curves for a single coverage metric on a given axis.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        Axis to plot on
    treatment_curves : dict
        Dictionary with treatment names as keys and tuples of (mean, min, max) curves as values
    title : str
        Title for the subplot
    colors : dict, optional
        Dictionary mapping treatment names to colors
    """
    if colors is None:
        colors = {
            'Unamended': '#7bccc4',
            'Catechin': '#fe9929'
        }
    
    max_samples = 0
    
    # Plot each treatment
    for treatment_name, (mean_curve, min_curve, max_curve) in treatment_curves.items():
        color = colors.get(treatment_name, None)
        
        x = np.arange(1, len(mean_curve) + 1)
        ax.plot(x, mean_curve, label=treatment_name, color=color, linewidth=2)
        ax.fill_between(x, min_curve, max_curve, alpha=0.3, color=color)
        
        max_samples = max(max_samples, len(mean_curve))
    
    # Axis labels
    ax.set_xlabel("Number of samples", fontsize=12)
    ax.set_ylabel("vOTUs detected (cumulative)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Set x-ticks based on the maximum number of samples
    xtick_positions = np.arange(1, max_samples + 1, 5)
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_positions, fontsize=10)
    ax.tick_params(axis='y', labelsize=10)
    
    # Start x-axis at 1
    ax.set_xlim(left=1)
    
    # Add legend
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(alpha=0.2)

def create_four_panel_species_accumulation(df_1_gene_mapped_per_10kb, df_prop10, df_1_read_mapped, df_5_reads_mapped, 
                                          treatment_cols_dict, permutations_count=100, 
                                          save_path=None):
    """
    Create four separate species accumulation plots in a column layout.
    
    Parameters:
    -----------
    df_1_gene_mapped_per_10kb : pd.DataFrame
        DataFrame with ≥1 gene mapped per 10kb
    df_prop10 : pd.DataFrame
        DataFrame with >10% of contig covered
    df_1_read_mapped : pd.DataFrame
        DataFrame with ≥1 read mapped
    df_5_reads_mapped : pd.DataFrame
        DataFrame with ≥5 reads mapped to structural phage gene
    treatment_cols_dict : dict
        Dictionary mapping treatment names to column lists
    permutations_count : int
        Number of permutations for accumulation curves
    save_path : str, optional
        Base path for saving figures (without extension)
    """
    # Create figure with four subplots in a column
    fig, axes = plt.subplots(4, 1, figsize=(6, 20))
    
    # Define dataset names and corresponding DataFrames
    datasets = [
        ('≥1 gene with ≥1 read mapped per 10kb', df_1_gene_mapped_per_10kb),
        ('>10% of contig covered', df_prop10),
        ('≥1 gene mapped', df_1_read_mapped),
        ('≥5 reads mapped to structural phage gene \n& ≥10% horizontal coverage', df_5_reads_mapped)
    ]
    
    all_treatment_curves = {}
    
    # Plot each dataset in a separate subplot
    for i, (title, df) in enumerate(datasets):
        print(f"Processing {title}...")
        
        # Prepare treatment data for this DataFrame
        treatment_dfs = prepare_treatment_data_for_single_df(df, treatment_cols_dict)
        
        # Calculate accumulation curves for each treatment
        treatment_curves = {}
        
        for treatment_name, df_treatment in treatment_dfs.items():
            print(f"  Calculating accumulation curves for {treatment_name}...")
            mean_curve, min_curve, max_curve = species_accumulation(df_treatment, permutations_count)
            treatment_curves[treatment_name] = (mean_curve, min_curve, max_curve)
        
        all_treatment_curves[title] = treatment_curves
        
        # Plot on the corresponding axis
        plot_single_species_accumulation(axes[i], treatment_curves, title)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figures if path provided
    if save_path:
        plt.savefig(f"{save_path}.pdf", dpi=300, bbox_inches='tight')
        plt.savefig(f"{save_path}.png", dpi=300, bbox_inches='tight')
    
    plt.show()
    
    return all_treatment_curves


# Define which columns belong to which treatments
treatment_cols_dict = {
    'Unamended': [
        "STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004",
        "STM_0716_E_M_E025", "STM_0716_E_M_E027", "STM_0716_E_M_E033",
        "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
        "STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E058"
    ],
    'Catechin': [
        "STM_0716_E_M_E059", "STM_0716_E_M_E060", "STM_0716_E_M_E062",
        "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E070",
        "STM_0716_E_M_E071", "STM_0716_E_M_E072", "STM_0716_E_M_E121",
        "STM_0716_E_M_E122", "STM_0716_E_M_E123", "STM_0716_E_M_E129",
        "STM_0716_E_M_E130", "STM_0716_E_M_E131"
    ]
}

# Create the four-panel plot
all_curves = create_four_panel_species_accumulation(
    df_1_gene_mapped_per_10kb=df_1_gene_mapped_per_10kb,
    df_prop10=df_prop10,
    df_1_read_mapped=df_1_read_mapped,
    df_5_reads_mapped=df_5_reads_mapped,
    treatment_cols_dict=treatment_cols_dict,
    permutations_count=100,
    save_path="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1X_specaccum_active_by_treatment_four_panel"
)


# # Investigate what alpha and beta diversity for ≥1 read mapped and ≥10% coverage looks like

# In[5]:


def plot_diversity_stats(df, outfile, figfile):

    # keep only the catechin and unamended columns
    active_df = df[["vOTU",
                              "STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004",
                              "STM_0716_E_M_E025", "STM_0716_E_M_E027", "STM_0716_E_M_E033",
                              "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
                              "STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E058",
                              "STM_0716_E_M_E059", "STM_0716_E_M_E060", "STM_0716_E_M_E062",
                              "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E070",
                              "STM_0716_E_M_E071", "STM_0716_E_M_E072", "STM_0716_E_M_E121",
                              "STM_0716_E_M_E122", "STM_0716_E_M_E123", "STM_0716_E_M_E129",
                              "STM_0716_E_M_E130", "STM_0716_E_M_E131"]]

    active_df = active_df.sort_values(by='vOTU')
    
    # for beta diversity
    active_df.to_csv(f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/{outfile}', sep='\t', index=False)
    
    active_df = active_df.set_index('vOTU')

    # drop rows that are all zeroes
    active_df = active_df.loc[active_df.sum(axis=1) > 0]
    print(len(active_df))
    
    replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')

    # Read metadata
    metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")
    
    # Simplify treatment labels
    metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
    metadata = metadata.loc[metadata['treatment'] != 'T']
    
    # Create treatment_timepoint column
    treatment_timepoint = metadata['treatment'] + "_" + metadata['timepoint'].astype(str)
    active_df.columns = treatment_timepoint.values
    
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
    
    active_df.columns = make_unique(active_df.columns)
    
    # Transpose and convert to numeric matrix
    data_mtx = active_df.T
    data_mtx = data_mtx.apply(pd.to_numeric, errors='coerce')

    metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('day')[1]))
    metadata['shannon_diversity'] = data_mtx.apply(shannon, axis=1).values

    # prop reads mapping to vOTUs

    gff_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt', sep='\t')
    gff_df = gff_df[['vOTU', 'gene']]
    gff_df.head()
    
    # Read the TSV file without header
    counts = pd.read_csv(
        "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED_no0s.tsv",
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
    counts = counts.merge(gff_df, on='gene', how='left').dropna()
    counts_long = counts.melt(id_vars=['gene', 'vOTU'], var_name='Sample', value_name='num_reads_mapped')
    counts_long = counts_long.groupby(['vOTU', 'Sample']).agg({'num_reads_mapped': 'sum'}).reset_index()
    counts_long = counts_long.merge(replicate_frame, on='Sample', how='left')
    counts_long = counts_long.loc[counts_long['treatment'] != 'CT']
    
    # keep only active vOTUs
    counts_long = counts_long.loc[counts_long['vOTU'].isin(list(active_df.index))]
    reads_mapped_per_sample = counts_long.groupby(['Sample', 'treatment', 'day', 'replicate']).agg({'num_reads_mapped': 'sum'}).reset_index()
    reads_mapped_per_sample = reads_mapped_per_sample.merge(metadata[['Sample', 'total reads (R1+R2)']], on='Sample', how='left')
    reads_mapped_per_sample['prop_mapped'] = reads_mapped_per_sample['num_reads_mapped'] / reads_mapped_per_sample['total reads (R1+R2)'] * 100
    reads_mapped_per_sample['day'] = reads_mapped_per_sample['day'].astype(int)
    metadata = metadata.merge(reads_mapped_per_sample[['Sample', 'prop_mapped']], on='Sample', how='left')

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
    plt.savefig(figfile, dpi=300)
    plt.show()


# In[6]:


plot_diversity_stats(df_1_read_mapped, 'active_vOTUs_1_read_mapped_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_1_read_mapped_alpha_diversity.pdf')


# In[7]:


plot_diversity_stats(df_1_gene_mapped_per_10kb, 'active_vOTUs_1_gene_per_10kb_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_1_gene_per_10kb_alpha_diversity.pdf')


# In[8]:


plot_diversity_stats(df_5_reads_mapped, 'active_vOTUs_5_reads_mapped_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_5_reads_mapped_alpha_diversity.pdf')


# In[9]:


plot_diversity_stats(df_prop10, 'active_vOTUs_prop10_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_prop10_alpha_diversity.pdf')


# # Plot rank abundance of all vOTUs, then show which get filtered out by each cutoff
# Needs 1) vOTU, 2) max abundance, and 3) categorized by detection level

# In[10]:


melted_df = df_1_read_mapped.melt(id_vars='vOTU', var_name='Sample', value_name='abundance')
melted_df


# In[11]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[12]:


melted_df = melted_df.merge(replicate_frame, on='Sample', how='left')


# In[13]:


melted_df


# In[14]:


# Your existing aggregation code
abundance_agg = melted_df.groupby(['vOTU', 'treatment']).agg({'abundance': 'max'})
abundance_agg = abundance_agg.reset_index()

# Pivot to get separate columns for each treatment
abundance_pivot = abundance_agg.pivot_table(
    index='vOTU', 
    columns='treatment', 
    values='abundance', 
    fill_value=0
).reset_index()

# Rank based on catechin treatment abundance
abundance_pivot = abundance_pivot.sort_values(by='unamended', ascending=False).sort_values(by='catechin', ascending=False).reset_index(drop=True)
abundance_pivot['rank'] = abundance_pivot.index + 1

# Calculate log10 abundance for both treatments
abundance_pivot['log10_catechin'] = abundance_pivot['catechin'].apply(lambda x: np.log10(x + 1))
abundance_pivot['log10_unamended'] = abundance_pivot['unamended'].apply(lambda x: np.log10(x + 1))


# In[15]:


abundance_pivot


# In[16]:


# Merge with provirus metadata
vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
provirus_df = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()
abundance_pivot = abundance_pivot.merge(provirus_df, on='vOTU', how='left')

# Sort by provirus status first, then by catechin abundance
abundance_pivot = abundance_pivot.sort_values(by=['is_provirus_aggregated', 'catechin'], ascending=False)
abundance_pivot['rank'] = np.arange(len(abundance_pivot))


# In[17]:


# add membership data
lytic_active = set(df_5_reads_mapped.vOTU.unique())
one_gene_per_10kb = set(df_1_gene_mapped_per_10kb.vOTU.unique())
prop10 = set(df_prop10.vOTU.unique())
one_read_mapped = set(df_1_read_mapped.vOTU.unique())


# In[28]:


def assign_color(vOTU):
    if vOTU in lytic_active:
        return 'lytic_active'
        
    if vOTU in prop10:
        return '>10% coverage'
        
    if vOTU in one_gene_per_10kb:
        return 'one_gene_per_10kb'
        
    if vOTU in one_read_mapped:
        return 'one_read_mapped'
        
    else:
        return 'grey'


# In[29]:


abundance_pivot['group'] = abundance_pivot['vOTU'].apply(assign_color)
abundance_pivot = abundance_pivot.loc[abundance_pivot['group'] != 'grey']
# Sort by provirus status first, then by catechin abundance
abundance_pivot = abundance_pivot.sort_values(by=['is_provirus_aggregated', 'catechin'], ascending=False)
abundance_pivot['rank'] = np.arange(len(abundance_pivot))


# In[31]:


# Find the boundary between provirus and non-provirus for vertical line
provirus_boundary = abundance_pivot['is_provirus_aggregated'].sum() - 0.5

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 6),
                               gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

# Rank-abundance plot with two points per vOTU
ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_catechin'], 
            c='#FC9B2D', 
            marker='o', s=20, alpha=0.3, label='Catechin')

ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_unamended'], 
            c='#7ACBC3', 
            marker='o', s=20, alpha=0.3, label='Unamended')

# Add vertical line to separate provirus from non-provirus
ax1.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8, 
            label='Provirus/Non-provirus boundary')

ax1.set_ylabel('log10(max trimmed mean abundance + 1)', size=14)
ax1.legend()

plt.tight_layout()
plt.show()


# In[32]:


# Find the boundary between provirus and non-provirus for vertical line
provirus_boundary = abundance_pivot['is_provirus_aggregated'].sum() - 0.5

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8),
                               gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

# Rank-abundance plot with two points per vOTU
ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_catechin'], 
            c='#FC9B2D', 
            marker='o', s=20, alpha=0.3, label='Catechin')
ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_unamended'], 
            c='#7ACBC3', 
            marker='o', s=20, alpha=0.3, label='Unamended')

# Add vertical line to separate provirus from non-provirus
ax1.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8, 
            label='Provirus/Non-provirus boundary')
ax1.set_ylabel('log10(max trimmed mean abundance + 1)', size=14)
ax1.legend()

# Heatmap bar chart using 'group' column for colors
# Create a colormap from the unique groups with specific colors: red, blue, black
unique_groups = abundance_pivot['group'].unique()
colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00']  # blue, orange, teal, vermillion
group_color_map = dict(zip(unique_groups, colors[:len(unique_groups)]))

# Create bar colors based on group values
bar_colors = [group_color_map[group] for group in abundance_pivot['group']]

# Create heatmap-style bar chart
bars = ax2.bar(abundance_pivot['rank'], 
               height=1,  # Fixed height for heatmap appearance
               width=1,   # Width of each bar
               color=bar_colors,
               edgecolor='none',
               alpha=0.8)

# Add vertical line to match the scatter plot
ax2.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8)

# Format the heatmap subplot
ax2.set_ylabel('Group', size=12)
ax2.set_xlabel('Rank', size=14)
ax2.set_ylim(0, 1)
ax2.set_yticks([0.5])
ax2.set_yticklabels(['Groups'])

# Create legend for groups
legend_elements = [plt.Rectangle((0,0),1,1, facecolor=group_color_map[group], 
                                label=group) for group in unique_groups]
ax2.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))

plt.tight_layout()
plt.show()


# In[33]:


# Find the boundary between provirus and non-provirus for vertical line
provirus_boundary = abundance_pivot['is_provirus_aggregated'].sum() - 0.5

# Create colormap for groups
unique_groups = abundance_pivot['group'].unique()
colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00']  # blue, orange, teal, vermillion
group_color_map = dict(zip(unique_groups, colors[:len(unique_groups)]))

# Create point colors based on group values
point_colors = [group_color_map[group] for group in abundance_pivot['group']]

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

# Catechin plot (top)
ax1.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_catechin'], 
            c=point_colors, 
            marker='o', s=20, alpha=0.7)

# Add vertical line to separate provirus from non-provirus
ax1.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8, 
            label='Provirus/Non-provirus boundary')
ax1.set_ylabel('log10(max trimmed mean abundance + 1)\n(Catechin)', size=14)
ax1.set_title('Catechin Treatment', size=16)
ax1.legend()

# Unamended plot (bottom)
ax2.scatter(abundance_pivot['rank'], 
            abundance_pivot['log10_unamended'], 
            c=point_colors, 
            marker='o', s=20, alpha=0.7)

# Add vertical line to separate provirus from non-provirus
ax2.axvline(x=provirus_boundary, color='gray', linestyle='--', alpha=0.8, 
            label='Provirus/Non-provirus boundary')
ax2.set_ylabel('log10(max trimmed mean abundance + 1)\n(Unamended)', size=14)
ax2.set_xlabel('Rank', size=14)
ax2.set_title('Unamended Control', size=16)
ax2.legend()

# Create shared legend for groups
legend_elements = [plt.Rectangle((0,0),1,1, facecolor=group_color_map[group], 
                                label=group) for group in unique_groups]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98), 
           title='Groups', title_fontsize=12)

plt.tight_layout()
plt.show()


# In[35]:


# Sort by group and then re-rank within each group by abundance
# First, let's create separate rankings for catechin and unamended within each group
abundance_pivot_sorted = abundance_pivot.sort_values(['group', 'log10_catechin'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted['catechin_rank'] = range(1, len(abundance_pivot_sorted) + 1)

abundance_pivot_sorted = abundance_pivot_sorted.sort_values(['group', 'log10_unamended'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted['unamended_rank'] = range(1, len(abundance_pivot_sorted) + 1)

# Find boundaries between groups
group_boundaries = []
# unique_groups = abundance_pivot_sorted['group'].unique()
unique_groups = ['lytic_active', '>10% coverage', 'one_gene_per_10kb', 'one_read_mapped']

cumsum = 0
for group in unique_groups:
    group_count = (abundance_pivot_sorted['group'] == group).sum()
    cumsum += group_count
    if cumsum < len(abundance_pivot_sorted):  # Don't add boundary after last group
        group_boundaries.append(cumsum + 0.5)

# Create colormap for groups
colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00']  # blue, orange, teal, vermillion
group_color_map = dict(zip(unique_groups, colors[:len(unique_groups)]))

# Create point colors based on group values
point_colors = [group_color_map[group] for group in abundance_pivot_sorted['group']]

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

# Catechin plot (top) - use catechin-based ranking
ax1.scatter(abundance_pivot_sorted['catechin_rank'], 
            abundance_pivot_sorted['log10_catechin'], 
            c=point_colors, 
            marker='o', s=20, alpha=0.7)

# Add vertical lines to separate groups
for boundary in group_boundaries:
    ax1.axvline(x=boundary, color='gray', linestyle='--', alpha=0.8)

ax1.set_ylabel('log10(max(GeTMM) + 1)', size=14)
ax1.set_title('Catechin', size=14)

# Unamended plot (bottom) - use unamended-based ranking
ax2.scatter(abundance_pivot_sorted['unamended_rank'], 
            abundance_pivot_sorted['log10_unamended'], 
            c=point_colors, 
            marker='o', s=20, alpha=0.7)

# Add vertical lines to separate groups
for boundary in group_boundaries:
    ax2.axvline(x=boundary, color='gray', linestyle='--', alpha=0.8)

ax2.set_ylabel('log10(max(geTMM) + 1)', size=14)
ax2.set_xlabel('Rank (within group)', size=14)
ax2.set_title('Unamended', size=14)

# Create shared legend for groups
legend_elements = [plt.Rectangle((0,0),1,1, facecolor=group_color_map[group], 
                                label=group) for group in unique_groups]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98), 
           title='Groups', title_fontsize=14).remove()

# Add group labels at the top
ax1_twin = ax1.twiny()
ax1_twin.set_xlim(ax1.get_xlim())

# Calculate midpoints of each group for labeling
group_midpoints = []
cumsum = 0
for group in unique_groups:
    group_count = (abundance_pivot_sorted['group'] == group).sum()
    midpoint = cumsum + group_count / 2
    group_midpoints.append(midpoint)
    cumsum += group_count

ax1_twin.set_xticks(group_midpoints)
ax1_twin.set_xticklabels(unique_groups, size=14)
# ax1_twin.set_xlabel('Group', size=14)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_rank_abundance.pdf', dpi=300)

plt.show()


# In[37]:


# Sort by group and then re-rank within each group by abundance
# First, let's create separate rankings for catechin and unamended within each group
abundance_pivot_sorted_catechin = abundance_pivot.sort_values(['group', 'log10_catechin'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted_catechin['catechin_rank'] = range(1, len(abundance_pivot_sorted_catechin) + 1)

abundance_pivot_sorted_unamended = abundance_pivot.sort_values(['group', 'log10_unamended'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted_unamended['unamended_rank'] = range(1, len(abundance_pivot_sorted_unamended) + 1)

# Define group order
unique_groups = ['lytic_active', '>10% coverage', 'one_gene_per_10kb', 'one_read_mapped']

# Debug: Check counts for each group
print("Catechin sorting - group counts:")
for group in unique_groups:
    count = (abundance_pivot_sorted_catechin['group'] == group).sum()
    print(f"  {group}: {count}")

print("\nUnamended sorting - group counts:")
for group in unique_groups:
    count = (abundance_pivot_sorted_unamended['group'] == group).sum()
    print(f"  {group}: {count}")

# Find boundaries between groups for each sorting
def calculate_group_boundaries(df, unique_groups):
    """Calculate boundaries between groups for a given dataframe"""
    boundaries = []
    cumsum = 0
    for group in unique_groups:
        group_count = (df['group'] == group).sum()
        cumsum += group_count
        if cumsum < len(df):  # Don't add boundary after last group
            boundaries.append(cumsum + 0.5)
    return boundaries, cumsum

group_boundaries_catechin, total_catechin = calculate_group_boundaries(abundance_pivot_sorted_catechin, unique_groups)
group_boundaries_unamended, total_unamended = calculate_group_boundaries(abundance_pivot_sorted_unamended, unique_groups)

print(f"\nCatechin boundaries: {group_boundaries_catechin}")
print(f"Unamended boundaries: {group_boundaries_unamended}")

# Create colormap for groups
colors = ['#0173B2', '#DE8F05', '#029E73', '#D55E00']  # blue, orange, teal, vermillion
group_color_map = dict(zip(unique_groups, colors[:len(unique_groups)]))

# Create point colors based on group values
point_colors_catechin = [group_color_map[group] for group in abundance_pivot_sorted_catechin['group']]
point_colors_unamended = [group_color_map[group] for group in abundance_pivot_sorted_unamended['group']]

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=False)

# Catechin plot (top) - use catechin-based ranking
ax1.scatter(abundance_pivot_sorted_catechin['catechin_rank'], 
            abundance_pivot_sorted_catechin['log10_catechin'], 
            c=point_colors_catechin, 
            marker='o', s=20, alpha=0.7)

# Add vertical lines to separate groups
for boundary in group_boundaries_catechin:
    ax1.axvline(x=boundary, color='gray', linestyle='--', alpha=0.8, linewidth=2)

ax1.set_ylabel('log10(max(GeTMM) + 1)', size=14)
ax1.set_title('Catechin', size=14)

# Unamended plot (bottom) - use unamended-based ranking
ax2.scatter(abundance_pivot_sorted_unamended['unamended_rank'], 
            abundance_pivot_sorted_unamended['log10_unamended'], 
            c=point_colors_unamended, 
            marker='o', s=20, alpha=0.7)

# Add vertical lines to separate groups
for boundary in group_boundaries_unamended:
    ax2.axvline(x=boundary, color='gray', linestyle='--', alpha=0.8, linewidth=2)

ax2.set_ylabel('log10(max(geTMM) + 1)', size=14)
ax2.set_xlabel('Rank (within group)', size=14)
ax2.set_title('Unamended', size=14)

# Create shared legend for groups
legend_elements = [plt.Rectangle((0,0),1,1, facecolor=group_color_map[group], 
                                label=group) for group in unique_groups]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98), 
           title='Groups', title_fontsize=14)

# Add group labels at the top (using catechin sorting)
ax1_twin = ax1.twiny()
ax1_twin.set_xlim(ax1.get_xlim())

# Calculate midpoints of each group for labeling
group_midpoints = []
cumsum = 0
for group in unique_groups:
    group_count = (abundance_pivot_sorted_catechin['group'] == group).sum()
    midpoint = cumsum + group_count / 2
    group_midpoints.append(midpoint)
    cumsum += group_count

ax1_twin.set_xticks(group_midpoints)
ax1_twin.set_xticklabels(unique_groups, size=14)

plt.tight_layout()
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3_rank_abundance.pdf', dpi=300)
plt.show()


# In[1]:


get_ipython().system('jupyter nbconvert --to script 002S3-active-virus-alternatives.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S3-AB_vOTU_specaccum_alpha_diversity_rank_abundance')


# In[ ]:




