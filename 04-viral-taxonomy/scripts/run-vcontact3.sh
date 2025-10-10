#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH --mem=180gb
#SBATCH -A PAS1117
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

vOTUFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna
inFile=../data/combined_manual_filtered_no_pipes.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-out


# remove pipe characters from input file because vcontact3 does not like them
cat ${vOTUFile} | sed -e 's/|/\_/g' > ${inFile}

mkdir -p ${outDir}

apptainer exec /fs/ess/PAS1117/modules/singularity/vcontact3.sif \
vcontact3 run \
--nucleotide ${inFile} \
--output ${outDir} \
--db-domain "prokaryotes" \
--db-version 228 \
--db-path /fs/ess/PAS1117/modules/sequence_dbs/vcontact3_dbs \
--exports cytoscape pyupset

# put the pipe characters back so outputs can be merged with other files (change '_provirus' back to '|provirus' to match vOTU file)
cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-out/exports/final_assignments.csv | sed -e 's/_provirus/|provirus/g' > /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-out/exports/final_assignments_with_pipe_chars.csv
cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-viral-taxonomy/results/vcontact3-out/exports/cytoscape/graph.cyjs | sed -e 's/_provirus/|provirus/g' > /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/04-taxonomy/results/vcontact3-out/exports/cytoscape/graph_with_pipe_chars.cyjs
