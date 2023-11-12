#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:30:30
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/noise_job.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/noise.%A_%a.err
#SBATCH --array=1-112%30 # Total tasks = 14 c_noise * 8 scale_noise

# Define arrays of values
c_noise_values=(0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85)
scale_noise_values=(1.0 2.0 2.5 3.0 3.5 4.3 4.5 5.0)

# Total number of c_noise values
num_c_noise=${#c_noise_values[@]}

# Calculate indices for c_noise and scale_noise
array_index=$(($SLURM_ARRAY_TASK_ID - 1))
c_noise_index=$(($array_index % num_c_noise))
scale_noise_index=$(($array_index / num_c_noise))

# Get the values for the current task
c_noise=${c_noise_values[$c_noise_index]}
scale_noise=${scale_noise_values[$scale_noise_index]}

# Run the Julia script with the current values
julia run_5.jl --T 4000 --c_noise $c_noise --scale_noise $scale_noise

