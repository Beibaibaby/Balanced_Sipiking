#!/bin/bash

# Define the main directory where all the results are stored
ARRAY_JOB_ID=$1  # Pass the ARRAY_JOB_ID as a command-line argument
main_dir="/gpfs/data/doiron-lab/draco/results_corr/exp_${ARRAY_JOB_ID}"

# Check if all tasks have completed
num_completed=$(ls "${main_dir}" | grep 'task_.*_completed' | wc -l)

if [ "$num_completed" -eq 80 ]; then
    # All tasks have completed
    julia plot_parallel.jl --main_dir "$main_dir"
else
    echo "Not all tasks have completed yet."
fi
