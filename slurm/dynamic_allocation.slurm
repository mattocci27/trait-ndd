#!/bin/bash
#SBATCH --job-name=dynamic_allocajtion
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=36:00:00
#SBATCH --array=ARRAY_INDICES_PLACEHOLDER
#SBATCH --output=log/model_%A_%a.out
#SBATCH --error=log/model_%A_%a.err

# ARRAY_INDICES="1-18" sbatch slurm/dynamic_allocation.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Pass the SLURM_ARRAY_TASK_ID to the R script to determine which model to run
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/dynamic_allocation.R ${SLURM_ARRAY_TASK_ID}

date
