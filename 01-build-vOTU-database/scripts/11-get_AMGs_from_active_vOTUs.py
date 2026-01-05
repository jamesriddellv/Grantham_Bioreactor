#!/usr/bin/env python
# coding: utf-8

# # The purpose of this notebook is to explore the AMGs identified from ALL vOTUs identified, then look at the subset that were expressed across clusters and determine significance. 

# We will try to interpret if any contribute to meaningful metabolisms, such as fermentation pathways, Nitrogen metabolism, and central carbon metabolism

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re


# In[2]:


dramv_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/DRAMv/DRAMv-annotate/annotations.tsv', sep='\t')
# assign column types to reduce memory
dramv_df = dramv_df.convert_dtypes()
dramv_df['gene_position'] = dramv_df['gene_position'].astype(int).astype(str)
dramv_df = dramv_df.rename(columns={'Unnamed: 0': 'gene'})
dramv_df['gene'] = dramv_df['gene'].str.replace(r'-cat_\d', '', regex=True)
dramv_df['gene'] = dramv_df['gene'].str.replace('_provirus', '|provirus')
dramv_df.head()


# ## Filter for active vOTUs

# In[3]:


dramv_df['vOTU'] = dramv_df['scaffold'].apply(lambda x: x.split('-cat')[0])
dramv_df['vOTU'] = dramv_df['vOTU'].apply(lambda x: x.replace('_provirus', '|provirus'))


# In[49]:


active_vOTUs = (
    list(
        pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/getmms_REV_5x_vOTUs_corrected_wide.tsv', sep='\t')
        ['vOTU']
    )
)
active_vOTUs


# In[50]:


len(active_vOTUs)


# In[51]:


dramv_df = dramv_df.loc[dramv_df['vOTU'].isin(active_vOTUs)]
dramv_df.vOTU.nunique()


# # Filter for AMGs (M flags) between V flags 

# In[52]:


import pandas as pd
import numpy as np

def find_votus_with_m_between_v(dramv_df):
    """
    Find vOTUs that contain a gene with 'M' flag (without 'F' or 'B') 
    positioned between two genes with 'V' flags.
    
    Parameters:
    dramv_df: DataFrame with columns ['vOTU', 'gene_position', 'amg_flags']
    
    Returns:
    tuple: (list of valid vOTUs, DataFrame of valid M genes)
    """
    
    def has_flag(flags_str, target_flag):
        """Check if a flag string contains a specific flag"""
        if pd.isna(flags_str):
            return False
        return target_flag in str(flags_str)
    
    def is_valid_m_gene(flags_str):
        """Check if gene has M flag but not F or B flags"""
        if pd.isna(flags_str):
            return False
        flags = str(flags_str)
        return 'M' in flags and 'F' not in flags and 'B' not in flags
    
    def has_v_flag(flags_str):
        """Check if gene has V flag"""
        if pd.isna(flags_str):
            return False
        return 'V' in str(flags_str)
    
    # Group by vOTU and process each group
    valid_votus = []
    valid_m_genes = []
    
    for votu, group in dramv_df.groupby('vOTU'):
        # Convert gene_position to integer and sort to ensure correct order
        group['gene_position_int'] = group['gene_position'].astype(int)
        group_sorted = group.sort_values('gene_position_int').reset_index(drop=True)
        
        # Get positions of genes with V flags
        v_positions = []
        for idx, row in group_sorted.iterrows():
            if has_v_flag(row['amg_flags']):
                v_positions.append(idx)
        
        # Get positions of genes with valid M flags (M without F or B)
        m_positions = []
        for idx, row in group_sorted.iterrows():
            if is_valid_m_gene(row['amg_flags']):
                m_positions.append(idx)
        
        # Check if any M gene is between two V genes
        found_valid_pattern = False
        
        for m_pos in m_positions:
            # Find V genes before and after this M gene
            v_before = [v for v in v_positions if v < m_pos]
            v_after = [v for v in v_positions if v > m_pos]
            
            if v_before and v_after:
                found_valid_pattern = True
                # Store this valid M gene
                valid_m_genes.append(group_sorted.iloc[m_pos])
        
        if found_valid_pattern:
            valid_votus.append(votu)
    
    # Convert valid M genes to DataFrame
    valid_m_df = pd.DataFrame(valid_m_genes) if valid_m_genes else pd.DataFrame()
    
    return valid_votus, valid_m_df

def analyze_votus_detailed(dramv_df):
    """
    Provide detailed analysis of vOTUs showing the flag patterns
    """
    
    def has_v_flag(flags_str):
        if pd.isna(flags_str):
            return False
        return 'V' in str(flags_str)
    
    def is_valid_m_gene(flags_str):
        if pd.isna(flags_str):
            return False
        flags = str(flags_str)
        return 'M' in flags and 'F' not in flags and 'B' not in flags
    
    results = []
    
    for votu, group in dramv_df.groupby('vOTU'):
        # Convert gene_position to integer and sort to ensure correct order
        group['gene_position_int'] = group['gene_position'].astype(int)
        group_sorted = group.sort_values('gene_position_int').reset_index(drop=True)
        
        # Create a summary of the flag pattern
        flag_pattern = []
        v_positions = []
        m_positions = []
        
        for idx, row in group_sorted.iterrows():
            flags = row['amg_flags']
            gene_pos = row['gene_position']
            
            if has_v_flag(flags):
                flag_pattern.append(f"V(pos_{gene_pos})")
                v_positions.append(idx)
            elif is_valid_m_gene(flags):
                flag_pattern.append(f"M(pos_{gene_pos})")
                m_positions.append(idx)
            elif pd.notna(flags):
                flag_pattern.append(f"{flags}(pos_{gene_pos})")
            else:
                flag_pattern.append(f"-(pos_{gene_pos})")
        
        # Check for valid pattern
        has_valid_pattern = False
        for m_pos in m_positions:
            v_before = [v for v in v_positions if v < m_pos]
            v_after = [v for v in v_positions if v > m_pos]
            if v_before and v_after:
                has_valid_pattern = True
                break
        
        results.append({
            'vOTU': votu,
            'has_valid_pattern': has_valid_pattern,
            'total_genes': len(group_sorted),
            'v_count': len(v_positions),
            'valid_m_count': len(m_positions),
            'flag_pattern': ' -> '.join(flag_pattern[:10])  # Show first 10 for readability
        })
    
    return pd.DataFrame(results)

# Example usage:
# Assuming your dataframe is called dramv_df
# valid_votus, valid_m_genes_df = find_votus_with_m_between_v(dramv_df)
# print(f"Found {len(valid_votus)} vOTUs with valid M-between-V patterns:")
# print(valid_votus)
# 
# print(f"\nFound {len(valid_m_genes_df)} valid AMG genes:")
# print(valid_m_genes_df[['vOTU', 'gene_position', 'amg_flags', 'kegg_hit', 'ko_id', 'pfam_hits', 'cazy_hits']].to_string())

# For detailed analysis:
# detailed_results = analyze_votus_detailed(dramv_df)
# valid_detailed = detailed_results[detailed_results['has_valid_pattern'] == True]
# print(valid_detailed)


# In[53]:


valid_votus = find_votus_with_m_between_v(dramv_df)
print(f"Found {len(valid_votus[0])} vOTUs with valid M-between-V patterns:")
valid_votus[1]


# In[54]:


valid_votus[1].columns


# In[55]:


valid_votus[1]


# In[56]:


detailed_results = analyze_votus_detailed(dramv_df)
valid_detailed = detailed_results[detailed_results['has_valid_pattern'] == True]


# In[57]:


valid_detailed


# In[58]:


dramv_df.loc[dramv_df['vOTU'] == 'contig_1167674']


# In[59]:


amg_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/DRAMv/DRAMv-distill/amg_summary.tsv', sep='\t')
amg_df['vOTU'] = amg_df['scaffold'].apply(lambda x: x.split('-cat')[0])
amg_df['vOTU'] = amg_df['vOTU'].apply(lambda x: x.replace('_provirus', '|provirus'))
amg_df


# In[60]:


amg_df['gene'] = amg_df['gene'].apply(lambda x: x.split('-cat_1')[0] + x.split('-cat_1')[1])


# In[61]:


amg_df = amg_df.loc[amg_df['gene'].isin(valid_votus[1]['gene'])]
amg_df


# In[62]:


iphop_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/data/vOTUs_MAGs_no_cutoff_for_cytoscape.tsv', sep='\t')
iphop_df = iphop_df.drop_duplicates(subset=['vOTU', 'highest_host_tax_rank'])
iphop_df = iphop_df.groupby(['vOTU']).agg({'highest_host_tax_rank': ';'.join}).reset_index()
iphop_df.head()


# In[63]:


amg_df = amg_df.merge(iphop_df, on='vOTU', how='left')


# In[64]:


amg_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/DRAMv/DRAMv-distill/amg_summary_QC_with_host_prediction.tsv', sep='\t', index=False)


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 06-AMG_EDA.ipynb --output')

