#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --account=PAS1117
#SBATCH --job-name=get_tar_contents
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1
#SBATCH --output=slurm-%x_%j.out


tar -tf ../data/metaT/trimmed_100M.tar.gz > metaT.conf
