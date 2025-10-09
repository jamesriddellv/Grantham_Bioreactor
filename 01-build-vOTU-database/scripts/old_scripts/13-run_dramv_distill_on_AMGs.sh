#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --mem=92gb
#SBATCH --account=PAS1117
#SBATCH --job-name=dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

module use /fs/project/PAS1117/modulefiles
module load DRAM

# For AMGs from cluster members

time DRAM-v.py distill \
-i /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/DRAMv/Permissive_AMG_annotations.tsv \
-o "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/DRAMv/AMG-dramv-distill"
