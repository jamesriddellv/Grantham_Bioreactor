#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass1
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=0-38%10

bash ./01-run_vs2_pass1.sh ${SLURM_ARRAY_TASK_ID}
