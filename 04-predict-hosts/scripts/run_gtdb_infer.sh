#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -n 40
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=gtdb_infer
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load GTDB-Tk

dataDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/MAGs"

# bacteria
# gtdbtk de_novo_wf --genome_dir ${dataDir}/ --bacteria --outgroup_taxon p__Patescibacteria --out_dir /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/Grantham_MAGs_GTDB-tk_results/ --cpus 40 --force --extension fa

# archaea
gtdbtk de_novo_wf --genome_dir ${dataDir}/ --archaea --outgroup_taxon p__Undinarchaeota --out_dir /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/Grantham_MAGs_GTDB-tk_results/ --cpus 40 --force --extension fa
