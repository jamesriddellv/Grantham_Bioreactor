#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=0-38%10

bash ./02-run_checkv.sh ${SLURM_ARRAY_TASK_ID}
