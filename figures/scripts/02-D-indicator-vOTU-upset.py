#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import plot


# In[2]:


indicators = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/indicator_vOTUs_per_cluster_multi.csv')
indicators = indicators.loc[indicators['p_value'] <= 0.01][['vOTU', 'significant_clusters']]
indicators.head()


# In[3]:


# Create a dictionary to store sets for each cluster
cluster_sets = {}

# Iterate through each row
for index, row in indicators.iterrows():
    votu = row['vOTU']
    clusters = row['significant_clusters']
    
    # Handle the clusters string - split by comma if multiple clusters
    if ',' in str(clusters):
        cluster_list = [c.strip() for c in str(clusters).split(',')]
    else:
        cluster_list = [str(clusters).strip()]
    
    # Add the vOTU to each cluster set
    for cluster in cluster_list:
        if cluster not in cluster_sets:
            cluster_sets[cluster] = set()
        cluster_sets[cluster].add(votu)

# Convert to the format you want (5 separate sets for clusters 1-5)
cluster_1 = cluster_sets.get('1', set())
cluster_2 = cluster_sets.get('2', set())
cluster_3 = cluster_sets.get('3', set())
cluster_4 = cluster_sets.get('4', set())
cluster_5 = cluster_sets.get('5', set())


# In[4]:


cluster_sets = [cluster_1, cluster_2, cluster_3, cluster_4, cluster_5]


# In[7]:


set_names = ['1', '2', '3', '4', '5']
all_elems = list(set().union(*cluster_sets))
df = pd.DataFrame([[e in st for st in cluster_sets] for e in all_elems], columns = set_names)
df_up = df.groupby(set_names).size()
plot(df_up, orientation='horizontal')
plt.savefig('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/02-D_indicator_upset.pdf', dpi=300)
plt.show()


# In[1]:


get_ipython().system('jupyter nbconvert --to script 005-indicator-vOTU-analysis.ipynb --output /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/scripts/02-D-indicator-vOTU-upset')


# In[ ]:




