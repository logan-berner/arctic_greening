#!/bin/bash
#SBATCH --job-name=calc_shrub_rwi
#SBATCH --output=/scratch/lb968/calc_shrub_rwi_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=00:02:00
#SBATCH --mem=750
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

echo calc shrub rwi
date

module load R/3.6.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/3.1.1_calc_shrub_rwi.R ${SLURM_ARRAY_TASK_ID} 
