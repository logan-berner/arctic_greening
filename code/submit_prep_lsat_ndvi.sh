#!/bin/bash
#SBATCH --job-name=prep_lsat
#SBATCH --output=/scratch/lb968/prep_lsat_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=03:00:00
#SBATCH --mem=20000
#SBATCH --cpus-per-task=1
#SBATCH --array=545,547,644-721,806-813,917,961

echo prep lsat ndvi
date

module load R/3.5.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/1.2_prep_lsat_ndvi_timeseries.R ${SLURM_ARRAY_TASK_ID} 
