#!/bin/bash

# Define an array with your desired SLURM array ranges
# declare -a ARRAY_RANGES=("1-18" "19-36" "37-48")
declare -a ARRAY_RANGES=("1-9" "10-18" "19-27" "28-36" "37-48")
# declare -a ARRAY_RANGES=("1-4")

# Loop through the array ranges and submit each one as a job
for RANGE in "${ARRAY_RANGES[@]}"
do
    echo "Submitting job for array indices $RANGE"

    # Use sed to replace the placeholder with the actual range
    sed "s/ARRAY_INDICES_PLACEHOLDER/${RANGE}%18/g" slurm/dynamic_allocation.slurm > temp.slurm

    # Submit the job using the temporary modified script
    sbatch temp.slurm

    if [ $? -eq 0 ]; then
        echo "Job submitted, sleeping for a while..."
        # Sleep for a specified duration, e.g., 5 minutes (300 seconds)
        sleep 300
    else
        echo "Failed to submit job for array indices $RANGE"
        # Exiting if sbatch fails
        exit 1
    fi
done

# Clean up the temporary script
rm temp.slurm

echo "All jobs submitted."
