#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:30:30
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/noisenof_job.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/noisenof.%A_%a.err
#SBATCH --array=1-150%50

# Define arrays of values
c_noise_values=(0.1 0.2 0.3 0.4 0.45 0.5 0.6 0.7 0.8 0.9)
scale_sigma_noises=(0.1 0.3 0.5 0.7 1.0 1.3 1.6 1.9 2.1 2.5 3.0 3.5 4.0 5.0 10.0)

# Total number of c_noise and scale_noise values
num_c_noise=${#c_noise_values[@]}
num_scale_noise=${#scale_sigma_noises[@]}

# Calculate indices for c_noise and scale_noise
array_index=$(($SLURM_ARRAY_TASK_ID - 1))
c_noise_index=$(($array_index % num_c_noise))
scale_noise_index=$(($array_index / num_c_noise))

# Ensure scale_noise_index is within the range of scale_sigma_noises
scale_noise_index=$(($scale_noise_index % num_scale_noise))

# Get the values for the current task
c_noise=${c_noise_values[$c_noise_index]}
scale_noise=${scale_sigma_noises[$scale_noise_index]}

# Run the Julia script with the current values
julia run_6.jl --T 2000 --c_noise $c_noise --sigma_noise $scale_noise


