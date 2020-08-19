#!/bin/bash
#SBATCH --job-name=calc_swi_ens
#SBATCH --output=/scratch/lb968/prep_swi_ensemble_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=00:35:00
#SBATCH --mem=400
#SBATCH --cpus-per-task=1
#SBATCH --array=1-35

echo prep  swi ensemble
date

module load R/3.5.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/2.2.1_prep_swi_ensemble.R ${SLURM_ARRAY_TASK_ID} 
