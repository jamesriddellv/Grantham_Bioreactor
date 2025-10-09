import pandas as pd
from Bio import SeqIO

def clean_fasta_headers(input_fasta, output_fasta):
    with open(output_fasta, "w") as out_f:
        headers_count = {}  # To track occurrences of each prefix

        # First pass: count occurrences of each prefix
        for record in SeqIO.parse(input_fasta, "fasta"):
            # remove everything after the first space
            without_checkv_mods = record.id.split(" ")[0]

            # extract the header before the last underscore that is added by checkV
            header_parts = without_checkv_mods.rsplit("_", 1)
            prefix = header_parts[0]
            headers_count[prefix] = headers_count.get(prefix, 0) + 1

        # Second pass: modify headers accordingly
        for record in SeqIO.parse(input_fasta, "fasta"):
            # remove everything after the first space
            without_checkv_mods = record.id.split(" ")[0]

            # extract the header before the last underscore that is added by checkV
            header_parts = without_checkv_mods.rsplit("_", 1)
            prefix = header_parts[0]

            if headers_count[prefix] == 1:
                new_header = prefix  # Remove last underscore and everything after

                record.id = new_header
                record.description = ""  # Remove description after space
                SeqIO.write(record, out_f, "fasta")

clean_fasta_headers("/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/proviruses.fna", "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/unmodified_proviruses.fna")
