#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH -n 15
#SBATCH --account=PAS1117
#SBATCH --job-name=download_iphop_db
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load iPHoP/1.3.3

iphop download --no_prompt --db_dir /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/iphop_Aug_2023_pub_rw
