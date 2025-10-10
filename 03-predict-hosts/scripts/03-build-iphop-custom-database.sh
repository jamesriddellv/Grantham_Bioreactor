#!/bin/bash
#SBATCH --time=00-12:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=iphop_add_to_db
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load iPHoP/1.3.3

iphop add_to_db --fna_dir /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/prefix/ \
--gtdb_dir ../data/prefix/Grantham_MAGs_GTDB-tk_results/ \
--out_dir ../data/prefix/iphop_Aug_2023_pub_prefix_grantham_rw \
--db_dir /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/iphop_Aug_2023_pub_rw/Aug_2023_pub_rw \
--num_threads 30
