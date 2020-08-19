#!/bin/bash
#SBATCH --job-name=calc_ndvi_trnds
#SBATCH --output=/scratch/lb968/calc_ndvi_trnds_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=00:12:00
#SBATCH --mem=3000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

echo calc lsat ndvi trends
date

module load R/3.5.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/2.1_calc_lsat_ndvi_trends.R ${SLURM_ARRAY_TASK_ID} 
