#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --account=PAS1117
#SBATCH --job-name=compress_vr-seq
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

tar cf - results | pigz -9 -p 20 > results.tar.gz
