import pandas as pd
import numpy as np
import os

data_dir = '../results/coverm-files'  # Replace with the path to your directory containing the files
print('reading files from ' + data_dir)

merged_df = []
for filename in os.listdir(data_dir):
    file_path = os.path.join(data_dir, filename)
    df = pd.read_csv(file_path, sep='\t')
    if type(merged_df) == list:
        merged_df = df
    else:
        merged_df = merged_df.merge(df, on='Contig', how='outer')

merged_df.to_csv(data_dir + '/ALL-SAMPLES-coverm_readcounts.txt', index=False)
print('merged coverage files to ' + data_dir + '/ALL-SAMPLES-coverm_readcounts.txt')
