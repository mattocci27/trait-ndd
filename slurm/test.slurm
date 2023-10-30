#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=log/test.out
#SBATCH --error=log/test.err
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=60:00:00

#Example:
# sbatch slurm/test.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Run the tar_make_clustermq() command with the specified number of clusters
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/run_script.R 3 > log/test.log

date