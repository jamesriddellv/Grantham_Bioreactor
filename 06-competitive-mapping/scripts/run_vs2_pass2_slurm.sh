#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

bash ./run_vs2_pass2.sh
