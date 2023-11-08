#!/bin/bash

# Define arrays of values
d_values=(0.05 0.1 0.3 0.5 0.7 1)
f_values=(0.0 0.2 0.5 0.7 0.9)
stimstr_values=(0.1 0.2 0.3 0.4 0.5 0.6)

# Make sure the trash directory exists
mkdir -p ../trash

# Define a function to run a single instance of the job with passed parameters
run_job() {
    local d=$1
    local f=$2
    local stimstr=$3
    local jobname="draco_test_${d}_${f}_${stimstr}"
    
    # Run the Julia command
    julia run_inject_3.0.0.jl --d $(printf "%.2f" $d) --f $(printf "%.2f" $f) --stimstr $(printf "%.2f" $stimstr) > "../trash/${jobname}.out" 2> "../trash/${jobname}.err" &
}

# Keep track of the number of concurrent jobs
max_jobs=20  # Set the max number of concurrent jobs
running_jobs=0

# Loop over parameter arrays
for d in "${d_values[@]}"; do
    for f in "${f_values[@]}"; do
        for stimstr in "${stimstr_values[@]}"; do
            # Run the job with the current set of parameters
            run_job $d $f $stimstr

            # Increment the count of running jobs
            ((running_jobs++))

            # Wait for jobs to finish if we've hit the max
            while [ $running_jobs -ge $max_jobs ]; do
                # Check the number of running jobs
                running_jobs=$(jobs -rp | wc -l)
                sleep 1
            done
        done
    done
done

# Wait for all jobs to finish
wait

echo "All experiments finished!"
