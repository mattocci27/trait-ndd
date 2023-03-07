#!/bin/bash
#SBATCH --job-name=model4
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=log/model4.out
#SBATCH --error=log/model4.err
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=96:00:00

#Example:
# sbatch slurm/model4.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Run the tar_make_clustermq() command with the specified number of clusters
# 18 jobs * 4 chains = 72 cpus
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/run_script.R 6 > log/model4.log

date
