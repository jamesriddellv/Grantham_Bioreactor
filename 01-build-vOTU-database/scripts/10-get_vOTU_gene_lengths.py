import re

# Define input and output file paths
input_gff = "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/prodigal-gv/combined_manual_filtered.gff"  # Change this to your actual file path
output_tsv = "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered_gene_lengths.txt"

with open(input_gff, "r") as infile, open(output_tsv, "w") as outfile:
    outfile.write("vOTU\tgene\tlength\n")  # Write header

    for line in infile:
        if line.startswith("#"):
            continue  # Skip comment lines

        fields = line.strip().split("\t")
        vOTU = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        length = end - start + 1  # Calculate gene length

        # Extract the gene ID from the ninth column
        match = re.search(r'ID=([^;]+)', fields[8])
        gene = match.group(1) if match else "NA"

        # Write to output file
        outfile.write(f"{vOTU}\t{gene}\t{length}\n")

print(f"Extraction complete! Results saved to {output_tsv}")
