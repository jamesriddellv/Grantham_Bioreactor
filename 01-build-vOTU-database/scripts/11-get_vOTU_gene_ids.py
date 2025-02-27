import re

# Input and output file paths
input_file = '../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/vOTU_genes.gff'
output_file = '../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/vOTU_gene_ids.tsv'

# Open the input file for reading
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Write the header to the output file
    outfile.write("vOTU\tgene\n")
    
    # Skip the first two lines
    lines = infile.readlines()[2:]
    
    # Process each remaining line
    for line in lines:
        # Skip lines that start with #
        if line.startswith('#'):
            continue
        
        parts = line.strip().split("\t")
        vOTU = parts[0]
        
        # Search for the gene ID using a regular expression
        gene_match = re.search(r'ID=([^;]+)', parts[-1])
        
        if gene_match:
            gene = gene_match.group(1)
            # Write the extracted values to the file
            outfile.write(f"{vOTU}\t{gene}\n")
        else:
            # Handle the case where the gene ID pattern wasn't found
            print(f"No gene ID found in line: {line}")

