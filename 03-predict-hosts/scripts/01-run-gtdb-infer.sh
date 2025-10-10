#!/bin/bash
#SBATCH --time=00-24:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=gtdb_infer
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

dataDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/prefix"

# bacteria
gtdbtk de_novo_wf \
--genome_dir ${dataDir}/ \
--bacteria \
--outgroup_taxon p__Patescibacteria \
--out_dir ../data/prefix_no_errors/Grantham_MAGs_GTDB-tk_results \
--cpus 28 \
--force \
--extension fa
wait

# archaea
gtdbtk de_novo_wf \
--genome_dir ${dataDir}/ \
--archaea \
--outgroup_taxon p__Undinarchaeota \
--out_dir ../data/prefix/Grantham_MAGs_GTDB-tk_results/ \
--cpus 28 \
--force \
--extension fa

cp gtdbtk.ar53.decorated.tree-table gtdbtk.ar122.decorated.tree-table
cp gtdbtk.ar53.decorated.tree gtdbtk.ar122.decorated.tree
