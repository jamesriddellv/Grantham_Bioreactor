import pandas as pd
import numpy as np

meancov_df = pd.read_csv('../results/ALL-SAMPLES-coverm_meancov.txt', sep=',')
cols = [name.strip(' Mean') for name in list(meancov_df.columns)]
meancov_df.columns = cols


total_reads = pd.read_csv('../results/coverm_stats.txt', sep=',')[['Sample', 'Total Mapped']]
total_reads['Sample'] = total_reads['Sample'].apply(lambda x: x.strip("\'"))


length = pd.read_csv('../results/vOTU_length.txt', sep='\t')
length.rename(columns={length.columns[1]: length.columns[1].split(' ')[1]}, inplace=True)

# sort everything then divide? instead of combining the dataframes.



meancov_length = meancov_df.merge(length, on='Contig', how='left')

length_normalized = pd.concat([meancov_length['Contig'], meancov_length.iloc[:,1:-1].div(meancov_length['Length'], axis=0)], axis=1)

# Assuming length_normalized and total_reads are your dataframes

# Set the "Contig" column as the index in length_normalized
length_normalized.set_index("Contig", inplace=True)

# Transpose the total_reads dataframe to make it easier to perform division
total_reads_transposed = total_reads.T
total_reads_transposed.columns = total_reads_transposed.iloc[0,:]
total_reads_transposed = total_reads_transposed.iloc[1,:]

# Divide length_normalized by the corresponding values in total_reads
result_df = length_normalized.divide(total_reads_transposed) * 100000000

# Reset the index if needed
result_df.reset_index(inplace=True)

result_df.to_csv('../results/vOTU_relative_abundance.tsv', sep='\t', index=False)
