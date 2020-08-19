#!/bin/bash
#SBATCH --job-name=calc_swi_trends
#SBATCH --output=/scratch/lb968/calc_swi_trends_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=00:18:00
#SBATCH --mem=3000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

echo calc swi ensemble
date

module load R/3.5.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/2.2.2_calc_swi_trends.R ${SLURM_ARRAY_TASK_ID} 
