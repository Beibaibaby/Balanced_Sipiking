#!/bin/bash

# Generate a unique identifier based on the current date and time
JOB_ID=$(date '+%Y%m%d_%H%M%S')

# Write the unique identifier to a file
echo $JOB_ID > unique_job_id

# Submit the SLURM job using the pl.sbatch script
sbatch pl.sbatch
