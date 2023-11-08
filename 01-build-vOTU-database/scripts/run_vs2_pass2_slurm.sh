#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=0,2-15

bash ./run_vs2_pass2.sh ${SLURM_ARRAY_TASK_ID}
