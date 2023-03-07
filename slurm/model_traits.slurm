#!/bin/bash
#SBATCH --job-name=model_traits
#SBATCH --partition=xhacnormalb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=log/model_traits.out
#SBATCH --error=log/model_traits.err
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=3000mb
#SBATCH --time=80:00:00

#Example:
# sbatch slurm/model_traits.slurm

date;hostname;pwd
module purge
module load singularity/3.7.3

# Run the tar_make_clustermq() command with the specified number of clusters
singularity exec --env RENV_PATHS_CACHE=${HOME}/renv --env RENV_PATHS_PREFIX_AUTO=TRUE radian.sif Rscript R/run_script.R $SLURM_CPUS_PER_TASK > log/model_traits.log

date
