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
#SBATCH --array=1-30%30 # 

# Define arrays of values
c_noise_values=(0.5 0.53 0.55 0.60 0.70 0.75)
scale_noise_values=(2.5 3.0 2.25 2.0 1.75)

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
julia run_5.jl --T 20000 --f_ee 0.0 --c_noise $c_noise --scale_noise $scale_noise

