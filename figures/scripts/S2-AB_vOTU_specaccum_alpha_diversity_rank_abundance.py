#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import shannon
from skbio.diversity.alpha import observed_otus
import seaborn as sns
import matplotlib as mpl
# set fonts and ensure PDF text is editable:
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'


from print_versions import print_versions
print_versions(globals())


# # This script investigates what species accumulation and alpha diversity would look like under different activity cutoffs

# In[2]:


df_1_read_mapped = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_vOTUs_wide.tsv', sep='\t')
df_prop10 = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1x_prop10_vOTUs_wide.tsv', sep='\t')
df_1_gene_mapped_per_10kb = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_1_gene_per_10kb_vOTUs_wide.tsv', sep='\t')
df_5_reads_mapped = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')


# # Species accumulation

# In[3]:


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
    Create four separate species accumulation plots in a column layout, ordered by stringency.
    """
    # Create figure with four subplots in a column
    fig, axes = plt.subplots(4, 1, figsize=(6, 20))
    
    # RE-ORDERED DATASETS LIST
    datasets = [
        ('≥1 read mapped', df_1_read_mapped),                                     # 1st: 1_read_mapped
        ('≥1 gene with ≥1 read mapped per 10kb', df_1_gene_mapped_per_10kb),      # 2nd: 1 gene with 1 read
        ('>10% of contig covered', df_prop10),                                    # 3rd: >10% coverage
        ('≥5 reads mapped to structural phage gene \n& ≥10% horizontal coverage', df_5_reads_mapped) # 4th: >5 reads
    ]
    
    all_treatment_curves = {}
    
    # Plot each dataset in the new order
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
    save_path="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S1-A_specaccum_active_by_treatment_four_panel"
)


# # Investigate what alpha and beta diversity for ≥1 read mapped and ≥10% coverage looks like

# In[4]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from skbio.diversity.alpha import hill

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
    
    # 1. ALIGNMENT STEP: Ensure metadata matches the EXACT order of active_df columns
    # Get the list of sample IDs currently in active_df
    sample_order = list(active_df.columns) 
    
    # Set Sample as index, reorder to match the dataframe, then reset index
    metadata = metadata.set_index('Sample').reindex(sample_order).reset_index()

    # 2. NOW simplify labels and calculate diversity
    metadata['treatment'] = metadata['treatment'].replace({'unamended': 'U', 'CT': 'T', 'catechin': 'C'})
    # Note: Ensure you don't drop 'T' samples if they are still in active_df, 
    # or drop them from BOTH active_df and metadata here.
    
    # Transpose and calculate diversity while names still match
    data_mtx = active_df.T.apply(pd.to_numeric, errors='coerce')
    metadata['hill_q1'] = data_mtx.apply(lambda x: hill(x, order=1), axis=1).values

    # 3. RENAME columns for plotting uniqueness AFTER diversity calculation
    treatment_timepoint = metadata['treatment'] + "_" + metadata['timepoint'].astype(str)
    
    def make_unique(columns):
        seen = {}
        new_cols = []
        for col in columns:
            seen[col] = seen.get(col, -1) + 1
            new_cols.append(col if seen[col] == 0 else f"{col}.{seen[col]}")
        return new_cols
    
    # Rename active_df columns only for your downstream mapping logic if needed
    active_df.columns = make_unique(treatment_timepoint.values)
    
    # Convert timepoint to int for plotting
    metadata['timepoint'] = metadata['timepoint'].apply(lambda x: int(x.split('day')[1]))

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
        "STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004", "STM_0716_E_M_E025",
        "STM_0716_E_M_E027", "STM_0716_E_M_E029", "STM_0716_E_M_E030", "STM_0716_E_M_E031",
        "STM_0716_E_M_E033", "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
        "STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E054", "STM_0716_E_M_E055",
        "STM_0716_E_M_E056", "STM_0716_E_M_E058", "STM_0716_E_M_E059", "STM_0716_E_M_E060",
        "STM_0716_E_M_E062", "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E066",
        "STM_0716_E_M_E067", "STM_0716_E_M_E068", "STM_0716_E_M_E070", "STM_0716_E_M_E071",
        "STM_0716_E_M_E072", "STM_0716_E_M_E121", "STM_0716_E_M_E122", "STM_0716_E_M_E123",
        "STM_0716_E_M_E125", "STM_0716_E_M_E126", "STM_0716_E_M_E127", "STM_0716_E_M_E129",
        "STM_0716_E_M_E130", "STM_0716_E_M_E131"
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

    # ── Z-Score Computations ──────────────────────────────────────────────────────
    metadata['hill_q1_z'] = metadata.groupby('timepoint')['hill_q1'].transform(
        lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
    )

    metadata['prop_mapped_log'] = np.log(metadata['prop_mapped'])
    metadata['prop_mapped_z'] = metadata.groupby('timepoint')['prop_mapped_log'].transform(
        lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
    )

    # ── Aesthetic mappings ────────────────────────────────────────────────────────
    color_palette = {"C": "#FC9B2D", "U": "#7ACBC3"}
    marker_map = {"0": "o", "7": "s", "14": "^", "21": "s", "35": "P"}

    # ── Helper functions for Stats and Plotting ───────────────────────────────────
    def get_groups(col):
        cat = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'C')][col].values
        una = metadata[(metadata['timepoint'] != 0) & (metadata['treatment'] == 'U')][col].values
        return cat, una

    def welch_one_tailed(col, h1_greater, label):
        """
        h1_greater : 'C' or 'U'  — which group is hypothesised to be greater
        """
        cat, una = get_groups(col)

        # Normality
        wc, pc = stats.shapiro(cat)
        wu, pu = stats.shapiro(una)
        print(f"\n=== Shapiro-Wilk [{label}] ===")
        print(f"  Catechin  W={wc:.4f} p={pc:.4f} {'✓' if pc > 0.05 else '✗'}")
        print(f"  Unamended W={wu:.4f} p={pu:.4f} {'✓' if pu > 0.05 else '✗'}")

        print(f"\n=== Descriptive Stats [{label}] ===")
        print(f"  Catechin  n={len(cat)} mean={np.mean(cat):.4f} SD={np.std(cat, ddof=1):.4f}")
        print(f"  Unamended n={len(una)} mean={np.mean(una):.4f} SD={np.std(una, ddof=1):.4f}")

        # One-tailed Welch's
        if h1_greater == 'U':
            t_stat, p_two = stats.ttest_ind(una, cat, equal_var=False)
            h1_label = "Unamended > Catechin"
        else:
            t_stat, p_two = stats.ttest_ind(cat, una, equal_var=False)
            h1_label = "Catechin > Unamended"

        p_one = p_two / 2 if t_stat > 0 else 1 - p_two / 2

        print(f"\n=== Welch's t-test (one-tailed: {h1_label}) [{label}] ===")
        print(f"  t = {t_stat:.4f}  p_two = {p_two:.4f}  p_one = {p_one:.4f}")
        sig = p_one < 0.05
        direction = h1_label if sig else f"No evidence that {h1_label}"
        print(f"  → {'Significant' if sig else 'Not significant'}: {direction}")

        return p_one

    def add_stripplot(ax, plot_df, y_col):
        for tp, m in marker_map.items():
            subset = plot_df[plot_df['timepoint_str'] == tp]
            if subset.empty:
                continue
            if tp == "7":
                sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                              marker='s', color='black', alpha=0.6, jitter=False, s=9)
                sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                              marker='x', color='white', alpha=1.0, jitter=False, s=6, linewidth=1)
            else:
                sns.stripplot(data=subset, x='treatment', y=y_col, ax=ax,
                              marker=m, color='black', alpha=0.6, jitter=True, s=8, linewidth=0.5)

    def get_stars(p):
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    def add_bracket(ax, y_vals, p_val):
        # Determine the symbol
        symbol = get_stars(p_val)
        
        # Calculate positioning
        y_max = np.nanmax(y_vals)
        y_pos = y_max + np.abs(y_max) * 0.05 # Used absolute value to handle potential negative z-scores
        h     = np.abs(y_max) * 0.03
        
        # Draw the bracket
        ax.plot([0, 0, 1, 1], [y_pos, y_pos+h, y_pos+h, y_pos], lw=1.2, c='k')
        
        # Add the text (stars or ns)
        # Shift 'ns' slightly higher or lower if needed for aesthetics
        v_offset = np.abs(y_max) * 0.01 if symbol == "ns" else 0 
        
        ax.text(0.5, y_pos + h + v_offset, symbol,
                ha='center', va='bottom', fontsize=12, 
                weight='bold' if symbol != "ns" else 'normal')

    # ── Run statistical tests on z-scores ─────────────────────────────────────────
    p_alpha = welch_one_tailed('hill_q1_z', h1_greater='U', label='Alpha Diversity (Z)')
    p_prop  = welch_one_tailed('prop_mapped_z', h1_greater='C', label='Prop Mapped (Z)')

    # ── Shared plot_df ────────────────────────────────────────────────────────────
    plot_df = metadata[~((metadata['treatment'] == 'U') & (metadata['timepoint'] == 0))].copy()
    plot_df['timepoint_str'] = plot_df['timepoint'].astype(str)

    # ── Build 2×2 figure ──────────────────────────────────────────────────────────
    sns.set_style("white")

    # Create a figure where the first column is 2x the width of the second
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), 
                             gridspec_kw={'width_ratios': [2, 0.5]})
    (ax_tl, ax_tr), (ax_bl, ax_br) = axes

    # ── TOP-LEFT: alpha diversity over time ───────────────────────────────────────
    sns.lineplot(data=metadata, x='timepoint', y='hill_q1_z',
                 hue='treatment', palette=color_palette,
                 estimator='mean', errorbar='ci', linewidth=1,
                 legend=True, ax=ax_tl)

    sns.scatterplot(data=metadata, x='timepoint', y='hill_q1_z',
                    hue='treatment', palette=color_palette,
                    legend=False, ax=ax_tl, s=50, alpha=0.6)

    ax_tl.set_xlabel("Day", fontsize=13)
    ax_tl.set_ylabel("vOTU Alpha Diversity Exp(H') Z-Score", fontsize=12)
    ax_tl.set_xticks([0, 7, 14, 21, 35])
    ax_tl.tick_params(axis='y', which='both', left=True, labelsize=11) 
    ax_tl.grid(True, alpha=0.3)
    for sp in ax_tl.spines.values(): sp.set_visible(False)

    # ── BOTTOM-LEFT: prop_mapped over time ────────────────────────────────────────
    sns.lineplot(data=metadata, x='timepoint', y='prop_mapped_z',
                 hue='treatment', palette=color_palette,
                 estimator='mean', errorbar='ci', linewidth=1, linestyle='--',
                 legend=True, ax=ax_bl)

    sns.scatterplot(data=metadata, x='timepoint', y='prop_mapped_z',
                    hue='treatment', palette=color_palette,
                    legend=False, ax=ax_bl, s=50, alpha=0.6)

    ax_bl.set_xlabel("Day", fontsize=13)
    ax_bl.set_ylabel("Log(% RNA mapped to active vOTUs) Z-Score", fontsize=12)
    ax_bl.set_xticks([0, 7, 14, 21, 35])
    ax_bl.tick_params(axis='y', which='both', left=True, labelsize=11)
    ax_bl.grid(True, alpha=0.3)
    for sp in ax_bl.spines.values(): sp.set_visible(False)

    # ── TOP-RIGHT: boxplot alpha diversity ────────────────────────────────────────
    sns.boxplot(data=plot_df, x='treatment', y='hill_q1_z',
                palette=color_palette, showfliers=False,
                width=0.6, linewidth=1.2, ax=ax_tr)
    add_stripplot(ax_tr, plot_df, 'hill_q1_z')
    add_bracket(ax_tr, plot_df['hill_q1_z'].values, p_alpha)
    ax_tr.set_xlabel("Treatment", fontsize=13)
    ax_tr.set_ylabel("vOTU Alpha Diversity Exp(H') Z-Score", fontsize=12)
    ax_tr.tick_params(axis='y', left=True)

    # ── BOTTOM-RIGHT: boxplot prop_mapped ─────────────────────────────────────────
    sns.boxplot(data=plot_df, x='treatment', y='prop_mapped_z',
                palette=color_palette, showfliers=False,
                width=0.6, linewidth=1.2, ax=ax_br)
    add_stripplot(ax_br, plot_df, 'prop_mapped_z')
    add_bracket(ax_br, plot_df['prop_mapped_z'].values, p_prop)
    ax_br.set_xlabel("Treatment", fontsize=13)
    ax_br.set_ylabel("Log(% RNA mapped to active vOTUs) Z-Score", fontsize=12)
    ax_br.tick_params(axis='y', left=True)

    # ── Save ──────────────────────────────────────────────────────────────────────
    fig.tight_layout()
    plt.savefig(figfile, dpi=300, bbox_inches='tight')
    plt.show()


# In[5]:


plot_diversity_stats(df_1_read_mapped, 'active_vOTUs_1_read_mapped_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2_1_read_mapped_alpha_diversity.pdf')


# In[6]:


plot_diversity_stats(df_1_gene_mapped_per_10kb, 'active_vOTUs_1_gene_per_10kb_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2_1_gene_per_10kb_alpha_diversity.pdf')


# In[7]:


plot_diversity_stats(df_5_reads_mapped, 'active_vOTUs_5_reads_mapped_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2_5_reads_mapped_alpha_diversity.pdf')


# In[8]:


plot_diversity_stats(df_prop10, 'active_vOTUs_prop10_relative_abundance.tsv', '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2_prop10_alpha_diversity.pdf')


# # Plot rank abundance of all vOTUs, then show which get filtered out by each cutoff
# Needs 1) vOTU, 2) max abundance, and 3) categorized by detection level

# In[9]:


melted_df = df_1_read_mapped.melt(id_vars='vOTU', var_name='Sample', value_name='abundance')
melted_df


# In[10]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')
replicate_frame.head()


# In[11]:


melted_df = melted_df.merge(replicate_frame, on='Sample', how='left')


# In[12]:


melted_df


# In[13]:


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


# In[14]:


abundance_pivot


# In[15]:


# Merge with provirus metadata
vOTU_metadata = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t')
provirus_df = vOTU_metadata[['vOTU', 'is_provirus_aggregated']].drop_duplicates()
abundance_pivot = abundance_pivot.merge(provirus_df, on='vOTU', how='left')

# Sort by provirus status first, then by catechin abundance
abundance_pivot = abundance_pivot.sort_values(by=['is_provirus_aggregated', 'catechin'], ascending=False)
abundance_pivot['rank'] = np.arange(len(abundance_pivot))


# In[16]:


# add membership data
lytic_active = set(df_5_reads_mapped.vOTU.unique())
one_gene_per_10kb = set(df_1_gene_mapped_per_10kb.vOTU.unique())
prop10 = set(df_prop10.vOTU.unique())
one_read_mapped = set(df_1_read_mapped.vOTU.unique())


# In[17]:


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


# In[18]:


abundance_pivot['group'] = abundance_pivot['vOTU'].apply(assign_color)
abundance_pivot = abundance_pivot.loc[abundance_pivot['group'] != 'grey']
# Sort by provirus status first, then by catechin abundance
abundance_pivot = abundance_pivot.sort_values(by=['is_provirus_aggregated', 'catechin'], ascending=False)
abundance_pivot['rank'] = np.arange(len(abundance_pivot))


# In[19]:


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


# In[20]:


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


# In[21]:


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


# In[22]:


# Define the desired order of groups
unique_groups = ['lytic_active', '>10% coverage', 'one_gene_per_10kb', 'one_read_mapped']

# Convert 'group' to categorical with specified order
abundance_pivot['group'] = pd.Categorical(abundance_pivot['group'], 
                                          categories=unique_groups, 
                                          ordered=True)

# Sort by group and then re-rank within each group by abundance
# First, let's create separate rankings for catechin and unamended within each group
abundance_pivot_sorted = abundance_pivot.sort_values(['group', 'log10_catechin'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted['catechin_rank'] = range(1, len(abundance_pivot_sorted) + 1)
abundance_pivot_sorted = abundance_pivot_sorted.sort_values(['group', 'log10_unamended'], ascending=[True, False]).reset_index(drop=True)
abundance_pivot_sorted['unamended_rank'] = range(1, len(abundance_pivot_sorted) + 1)

# Find boundaries between groups
group_boundaries = []
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

plt.tight_layout()
# plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2_rank_abundance.pdf', dpi=300)
plt.show()


# In[23]:


abundance_pivot_sorted


# In[24]:


abundance_pivot_sorted['group'].value_counts()


# In[25]:


annotations_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t', usecols=['vOTU', 'length', 'concat_annotations'])
annotations_df.head()


# In[26]:


annotations_grouped = annotations_df.groupby('vOTU').agg({'concat_annotations': ';'.join}).reset_index()


# In[27]:


abundance_pivot_sorted = abundance_pivot_sorted.merge(annotations_grouped, on='vOTU', how='left')


# In[28]:


lysis_keywords = [
    "endolysin", "holin", "antiholin", "spanin", "i-spanin", "o-spanin",
    "murein hydrolase", "peptidoglycan hydrolase",
    "N-acetylmuramoyl-L-alanine amidase", "lysozyme",
    "lytic transglycosylase", "muramidase", "autolysin", "pinholin",
    "Rz"
]

replication_keywords = [
    "DNA polymerase", "helicase", "primase", "DNA ligase",
    "single-stranded binding protein", "SSB", "replication protein",
    "rep protein", "DnaB", "DnaC", "DnaG", "origin binding protein",
    "OBP", "replicase", "DNA-dependent DNA polymerase", "sliding clamp",
    "clamp loader", "topoisomerase", "gyrase",
    "ribonuclease H", "RNaseH", "recombinase", "strand transferase",
    "exonuclease", "endonuclease"
]

integration_keywords = [
    # Integration / excision
    "integrase",
    "recombinase", "excisionase",
    "recombination directionality factor", "RDF", "phage integrase",
    "lysogeny", "repressor", "prophage", "ParB", "partition"
]

structural_keywords = ['capsid', 'tail', 'packag', 'portal', 'terminase', 'spike', 'baseplate', 'internal virion', 'neck', 'prohead']


def _annotation_matches(annotation_string, keywords):
    """
    Split a semicolon-delimited annotation string into individual annotations
    and check whether any annotation contains a keyword (case-insensitive).
    Returns True if any keyword is found in any annotation field.
    """
    if not isinstance(annotation_string, str) or annotation_string.strip() == "":
        return False
    annotations = annotation_string.split(";")
    for annotation in annotations:
        annotation_lower = annotation.lower()
        for keyword in keywords:
            if keyword.lower() in annotation_lower:
                return True
    return False


abundance_pivot_sorted["lysis"] = abundance_pivot_sorted["concat_annotations"].apply(
    lambda x: _annotation_matches(x, lysis_keywords)
)

abundance_pivot_sorted["replication"] = abundance_pivot_sorted["concat_annotations"].apply(
    lambda x: _annotation_matches(x, replication_keywords)
)

abundance_pivot_sorted["integration"] = abundance_pivot_sorted["concat_annotations"].apply(
    lambda x: _annotation_matches(x, integration_keywords)
)

abundance_pivot_sorted["structural"] = abundance_pivot_sorted["concat_annotations"].apply(
    lambda x: _annotation_matches(x, structural_keywords)
)



# In[29]:


abundance_pivot_sorted.columns


# In[30]:


abundance_pivot_sorted


# In[31]:


abundance_pivot_sorted[['lysis', 'replication', 'integration', 'structural']].sum()


# In[32]:


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 2. Prepare Heatmap Data
# Extract the boolean columns and transpose so categories are rows
anno_cols = ['lysis', 'replication', 'integration', 'structural']
# Convert True/False to 1/0 for plotting
heatmap_data = abundance_pivot_sorted[anno_cols].astype(int).T.values

# 3. Plotting
# Increase height_ratios so the heatmap is much shorter than the scatter plots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 12), sharex=True, 
                                    gridspec_kw={'height_ratios': [4, 4, 1.2]})

# --- Catechin Plot ---
ax1.scatter(abundance_pivot_sorted['catechin_rank'], 
            abundance_pivot_sorted['log10_catechin'], 
            c=point_colors, marker='o', s=20, alpha=0.7)
ax1.set_ylabel('log10(GeTMM + 1)', size=12)
ax1.set_title('Catechin Abundance', size=14, loc='left')

# --- Unamended Plot ---
ax2.scatter(abundance_pivot_sorted['catechin_rank'], 
            abundance_pivot_sorted['log10_unamended'], 
            c=point_colors, marker='o', s=20, alpha=0.7)
ax2.set_ylabel('log10(GeTMM + 1)', size=12)
ax2.set_title('Unamended Abundance', size=14, loc='left')

# --- Annotation Heatmap ---
# Added origin='lower' to put the first row at the bottom
im = ax3.imshow(heatmap_data, aspect='auto', cmap='Greys', interpolation='nearest',
                origin='lower',
                extent=[0.5, len(abundance_pivot_sorted) + 0.5, -0.5, len(anno_cols) - 0.5])

# Formatting remains the same
ax3.set_yticks(range(len(anno_cols)))
ax3.set_yticklabels(anno_cols, size=10)
# Formatting the Heatmap axis
ax3.set_yticks(range(len(anno_cols)))
ax3.set_yticklabels(anno_cols, size=10)
ax3.set_xlabel('vOTU Rank (ordered by group and Catechin abundance)', size=12)
ax3.tick_params(axis='x', which='both', bottom=False, labelbottom=False) # Hide x ticks for heatmap

# Add vertical group boundaries to ALL plots
for ax in [ax1, ax2, ax3]:
    for boundary in group_boundaries:
        ax.axvline(x=boundary, color='gray', linestyle='--', alpha=0.5, linewidth=1)

# Add group labels at the very top using a twin axis on ax1
ax1_twin = ax1.twiny()
ax1_twin.set_xlim(ax1.get_xlim())
ax1_twin.set_xticks(group_midpoints)
ax1_twin.set_xticklabels(unique_groups, size=12, fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(hspace=0.2) # Adjust spacing between subplots
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S2-B_rank_abundance.pdf', dpi=300)
plt.show()


# In[33]:


from venny4py.venny4py import *
import matplotlib.pyplot as plt

# Define the categories of interest
anno_cols = ['structural', 'lysis', 'integration', 'replication']

for group_name in abundance_pivot_sorted['group'].unique():
    
    group_df = abundance_pivot_sorted[abundance_pivot_sorted['group'] == group_name]
    
    # Map column names to sets of indices
    sets = {
        col: set(group_df[group_df[col] == True].index) 
        for col in anno_cols
    }
    
    if all(len(s) == 0 for s in sets.values()):
        continue

    # FIX: Removed 'as_no_gui'. Most versions only require 'sets' and 'out'.
    # If you want to show it in a notebook/IDE, set out='save' (saves to file)
    # or omit 'out' to just manipulate the plot directly.
    venny4py(sets=sets)
    
    plt.title(f"Functional Overlap: {group_name}", fontsize=14, pad=25)
    plt.show()


# # What proportion of abundance do structural-only contigs make up?

# In[34]:


structural_lytic_active_only = abundance_pivot_sorted.loc[
((abundance_pivot_sorted['structural'] == True)
 | (abundance_pivot_sorted['lysis'] == True)
)
& (abundance_pivot_sorted['integration'] == False)
& (abundance_pivot_sorted['replication'] == False)

& (abundance_pivot_sorted['group'] == 'lytic_active')
]


# In[35]:


lytic_active_sum_cat = abundance_pivot_sorted.loc[abundance_pivot_sorted['group'] == 'lytic_active']['catechin'].sum()
lytic_active_sum_unamended = abundance_pivot_sorted.loc[abundance_pivot_sorted['group'] == 'lytic_active']['unamended'].sum()


# In[36]:


structural_active_sum_cat = structural_lytic_active_only['catechin'].sum()
structural_active_sum_unamended = structural_lytic_active_only['unamended'].sum()


# In[37]:


structural_active_sum_cat / lytic_active_sum_cat


# In[38]:


structural_active_sum_unamended / lytic_active_sum_unamended


# In[39]:


# get length distributions
structural_lytic_active_only


# In[40]:


structural_vOTUs = structural_lytic_active_only.vOTU.unique().tolist()


# In[41]:


length_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', usecols=['vOTU', 'length'], sep='\t').drop_duplicates()
length_df


# In[42]:


sns.histplot(length_df.loc[length_df['vOTU'].isin(structural_vOTUs)], x='length', bins=30)


# # What are their host predictions

# In[43]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.sort_values(by=['vOTU', 'Confidence score']).drop_duplicates('vOTU')


# In[51]:


structural_iphop_df = iphop_df.loc[iphop_df['vOTU'].isin(structural_vOTUs)]


# In[54]:


structural_iphop_df.highest_host_tax_rank.value_counts()


# In[14]:


replicate_frame = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/sample_metadata.csv')

# Read metadata
metadata = pd.read_csv("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/metaT_sample_metadata.csv")


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
"STM_0716_E_M_E002", "STM_0716_E_M_E003", "STM_0716_E_M_E004", "STM_0716_E_M_E025",
"STM_0716_E_M_E027", "STM_0716_E_M_E029", "STM_0716_E_M_E030", "STM_0716_E_M_E031",
"STM_0716_E_M_E033", "STM_0716_E_M_E034", "STM_0716_E_M_E035", "STM_0716_E_M_E050",
"STM_0716_E_M_E051", "STM_0716_E_M_E052", "STM_0716_E_M_E054", "STM_0716_E_M_E055",
"STM_0716_E_M_E056", "STM_0716_E_M_E058", "STM_0716_E_M_E059", "STM_0716_E_M_E060",
"STM_0716_E_M_E062", "STM_0716_E_M_E063", "STM_0716_E_M_E064", "STM_0716_E_M_E066",
"STM_0716_E_M_E067", "STM_0716_E_M_E068", "STM_0716_E_M_E070", "STM_0716_E_M_E071",
"STM_0716_E_M_E072", "STM_0716_E_M_E121", "STM_0716_E_M_E122", "STM_0716_E_M_E123",
"STM_0716_E_M_E125", "STM_0716_E_M_E126", "STM_0716_E_M_E127", "STM_0716_E_M_E129",
"STM_0716_E_M_E130", "STM_0716_E_M_E131"
]
counts = counts.merge(gff_df, on='gene', how='left').dropna()
counts_long = counts.melt(id_vars=['gene', 'vOTU'], var_name='Sample', value_name='num_reads_mapped')
counts_long = counts_long.groupby(['vOTU', 'Sample']).agg({'num_reads_mapped': 'sum'}).reset_index()
counts_long = counts_long.merge(replicate_frame, on='Sample', how='left')
counts_long = counts_long.loc[counts_long['treatment'] != 'CT']

reads_mapped_per_sample = counts_long.groupby(['Sample', 'treatment', 'day', 'replicate']).agg({'num_reads_mapped': 'sum'}).reset_index()
reads_mapped_per_sample = reads_mapped_per_sample.merge(metadata[['Sample', 'total reads (R1+R2)']], on='Sample', how='left')
reads_mapped_per_sample['prop_mapped'] = reads_mapped_per_sample['num_reads_mapped'] / reads_mapped_per_sample['total reads (R1+R2)'] * 100
reads_mapped_per_sample['day'] = reads_mapped_per_sample['day'].astype(int)
metadata = metadata.merge(reads_mapped_per_sample[['Sample', 'prop_mapped']], on='Sample', how='left').dropna()


# In[20]:


np.std(metadata.prop_mapped)


# In[21]:


np.mean(metadata.prop_mapped)


# In[23]:


# 1. Define your probability thresholds (0.0, 0.1, ..., 1.0)
quantiles = np.linspace(0, 1, 11)

# 2. Calculate the quantiles
# If you want the quantiles of the raw values in that column:
deciles = np.quantile(counts_long['num_reads_mapped'], quantiles)

# Display results as a Series for better readability
quantile_table = pd.Series(deciles, index=[f"{int(q*100)}%" for q in quantiles])
print(quantile_table)


# In[45]:


get_ipython().system('jupyter nbconvert --to script 002S2-active-virus-alternatives.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/S2-AB_vOTU_specaccum_alpha_diversity_rank_abundance')


# In[ ]:





# In[ ]:




