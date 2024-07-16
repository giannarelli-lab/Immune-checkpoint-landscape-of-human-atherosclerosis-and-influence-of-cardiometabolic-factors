#!/bin/bash -l

#SBATCH -J cellchatr_conditions
#SBATCH --output="cellchatr_conditions-%A_%a.out"
#SBATCH -p cpu_medium
#SBATCH --export=ALL
#SBATCH --time=24:00:00
#SBATCH --mem=100000

module load r/4.0.3
Rscript cellchat.R