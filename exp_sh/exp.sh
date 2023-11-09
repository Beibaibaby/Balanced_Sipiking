#!/bin/bash

# Array of 'd' values
d_values=(0.05 0.10)

# Array of 'f' values
f_values=(0.05 0.1 0.2 0.3)

# Keep track of the number of cores
num_cores=$(nproc)
running_jobs=0

for d in "${d_values[@]}"
do
  for f in "${f_values[@]}"
  do
    echo "Running experiment with d = $d and f = $f on background"
    julia run_inject_3.0.0.jl --d $d --f $f --T 1000 --stim_start_time 500 --stim_duration 3 --stimstr 1.0 &

    # Increment the count of running jobs
    ((running_jobs++))

    # If we've reached the number of cores, wait for all jobs to complete
    if (( running_jobs == num_cores )); then
      wait # Wait for all background jobs to finish
      running_jobs=0 # Reset the running jobs counter
      echo "waiting"
    fi
  done
done

# Wait for the remaining jobs to finish, if any
wait

echo "All experiments finished!"