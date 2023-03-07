#!/bin/bash
#SBATCH --job-name=model3
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=log/model3.out
#SBATCH --error=log/model3.err
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=96:00:00

#Example:
# sbatch slurm/model3.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Run the tar_make_clustermq() command with the specified number of clusters
# 18 jobs * 4 chains = 72 cpus
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/run_script.R 18 > log/model3.log

date
