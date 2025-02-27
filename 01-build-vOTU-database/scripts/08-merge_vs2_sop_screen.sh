#!/bin/bash

workDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-sop-screen"

# merge scaffold outputs into one fasta, store in the permanent data directory
cat ${workDir}/*/final-viral-scored.fa > ${workDir}/vs2sop-final-viral-scored.fa
