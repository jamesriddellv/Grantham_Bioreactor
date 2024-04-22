#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -n 20
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=iphop_add_to_db
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load iPHoP/1.3.2

iphop add_to_db --fna_dir /fs/scratch/SullivanLab/JamesR/Grantham_Bioreactor/MAGs/ \
 --gtdb_dir /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/Grantham_MAGs_GTDB-tk_results/ \
--out_dir iphop_2021_pub_grantham_target_rw --db_dir /fs/project/PAS1117/modules/sequence_dbs/iPHoP/Sept_2021_pub_rw/
