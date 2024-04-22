import csv

def extract_cds_lengths(gff_file, output_file):
    # Open the GFF file for reading
    with open(gff_file, 'r') as file:
        # Initialize a CSV writer to write the output to the specified file
        with open(output_file, 'w', newline='') as csvfile:
            # Use a CSV writer with tab as the delimiter
            writer = csv.writer(csvfile, delimiter='\t')
            # Write the header row
            writer.writerow(["gene", "length"])
            
            # Process each line in the GFF file
            for line in file:
                # Skip comments and blank lines
                if line.startswith("#") or not line.strip():
                    continue
                
                # Split the line by tabs to get the fields
                fields = line.strip().split('\t')
                
                # Unpack the fields (assuming GFF format)
                seqname, source, feature, start, end, score, strand, frame, attribute = fields
                
                # Process only "CDS" features
                if feature == "CDS":
                    # Calculate the length of the CDS
                    length = int(end) - int(start) + 1
                    
                    # Extract gene ID from attributes field
                    # Assuming the gene ID is in the format "ID=gene_id" in the attributes field
                    # Modify this part to match your specific GFF file's attribute format if needed
                    attributes = attribute.split(';')
                    gene_id = None
                    for attr in attributes:
                        if attr.startswith("ID="):
                            gene_id = attr.split('=')[1]
                            break
                    
                    # If gene ID was found, write the gene ID and length to the output file
                    if gene_id:
                        writer.writerow([gene_id, length])
    
    print(f"CDS lengths extracted to: {output_file}")

# Specify the GFF file path and the output file path
gff_file = '../results/vOTU_clusters/checkv_5kb_or_high_quality_vOTUs.gff'  # Replace with the path to your GFF file
output_file = '../results/vOTU_clusters/5kb_or_high_quality_vOTU_gene_lengths.tsv'  # Specify the desired output file path (tab-separated file)

# Run the function to extract CDS lengths
extract_cds_lengths(gff_file, output_file)
