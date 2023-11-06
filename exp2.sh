#!/bin/bash

# Array of 'd' values
lambda_noises=(4.0 6.0 8.0 15.0)


# Keep track of the number of cores
num_cores=$(nproc)
running_jobs=0

for lambda_noise in "${lambda_noises[@]}"
do
  
    julia run_inject_3.0.0.jl --d 0.1 --f 0.5 --T 1000 --lambda_noise $lambda_noise &

    # Increment the count of running jobs
    ((running_jobs++))

    # If we've reached the number of cores, wait for all jobs to complete
    if (( running_jobs == num_cores )); then
      wait # Wait for all background jobs to finish
      running_jobs=0 # Reset the running jobs counter
      echo "waiting"
    fi
done

# Wait for the remaining jobs to finish, if any
wait

echo "All experiments finished!"