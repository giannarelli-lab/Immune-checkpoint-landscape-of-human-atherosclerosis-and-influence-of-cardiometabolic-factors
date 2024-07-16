#!/bin/bash -l

#SBATCH -J cellchat_conditions
#SBATCH --output="slurm/cellchat_conditions-%A_%a.out"
#SBATCH --array=1-12
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=8:00:00
#SBATCH --mem=150000
#SBATCH --exclusive

module load r/4.0.3

command=$(sed -n "$SLURM_ARRAY_TASK_ID"p commands.txt)
echo $command
srun $command