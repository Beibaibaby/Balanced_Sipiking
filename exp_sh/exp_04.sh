#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/jobname.%A_%a.err
#SBATCH --array=1-5%5

# Create the trash directory if it doesn't exist

# Define arrays of values
d_values=(0.3)
f_values=(0.9)
stimstr_values=(0.1 0.2 0.3 0.4 0.5)


# Total number of variations for each parameter
num_d=${#d_values[@]}
num_f=${#f_values[@]}
num_stimstr=${#stimstr_values[@]}

# Calculate the total number of combinations for d and f
total_df=$((num_d * num_f))

# Calculate parameter index for d, f, and stimstr
# ... rest of your script logic to determine the appropriate parameters ...

# Explicitly cast values as Float64 when passing to Julia script
julia run_inject_3.0.0.jl --d $(printf "%.2f" $d) --f $(printf "%.2f" $f) --stimstr $(printf "%.2f" $stimstr)
