#!/bin/bash

# Loop from 0 to 1 in increments of 0.1
for f in $(seq 0 0.1 1); do
    # Call your Julia script with different values of f
    julia run_inject_2.0.0.jl --T 10000 --f $f
done

# Print when all experiments are done
echo "All experiments finished!"