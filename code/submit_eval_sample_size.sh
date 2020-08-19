#!/bin/bash
#SBATCH --job-name=eval_sample_size
#SBATCH --output=/scratch/lb968/eval_sample_size_%a.log
#SBATCH --chdir=/scratch/lb968/
#SBATCH --time=01:10:00
#SBATCH --mem=2000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

echo eval lsat sample size
date

module load R/3.5.2
R --version

# run application
Rscript /home/lb968/code/arctic_greening/2.1_eval_lsat_ndvi_trends_sample_size.R ${SLURM_ARRAY_TASK_ID} 
