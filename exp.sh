#!/bin/bash

# Array of 'd' values
declare -a d_values=(0.01 0.10)

# Array of 'f' values
declare -a f_values=(0.05 0.1 0.2 0.3 0.40)

# Loop over 'd' values
for d in "${d_values[@]}"
do
  # Loop over 'f' values
  for f in "${f_values[@]}"
  do
    echo "Running experiment with d = $d and f = $f"
    julia run_inject_3.0.0.jl --d $d --f $f
  done
done

echo "All experiments finished!"