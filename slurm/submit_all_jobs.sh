#!/bin/bash

# Define an array with your desired SLURM array ranges
declare -a ARRAY_RANGES=("1-18" "19-36" "37-48")

# Loop through the array ranges and submit each one as a job
for RANGE in "${ARRAY_RANGES[@]}"
do
    echo "Submitting job for array indices $RANGE"
    ARRAY_INDICES="$RANGE" sbatch slurm/dynamic_allocation.slurm
    echo "Job submitted, sleeping for a while..."

    # Sleep for a specified duration, e.g., 5 minutes (300 seconds)
    sleep 10
done

echo "All jobs submitted."
