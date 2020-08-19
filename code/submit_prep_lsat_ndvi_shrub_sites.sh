#!/bin/bash
#SBATCH --job-name=prep_shrub_lsat
#SBATCH --output=/scratch/lb968/prep_shrub_lsat_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=2:30
#SBATCH --mem=1500
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

echo prep shrub lsat
date

module load R/3.6.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/3.1.2_prep_lsat_ndvi_for_shrub_rwi_comparison.R ${SLURM_ARRAY_TASK_ID} 
