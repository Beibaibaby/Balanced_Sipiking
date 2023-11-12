#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/noise_job.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/noise.%A_%a.err
#SBATCH --array=1-21%21

# Define an array of c_noise values (21 values from 0.0 to 1.0 in steps of 0.05)
c_noise_values=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)

# Get the c_noise value for the current task
c_noise_index=$(($SLURM_ARRAY_TASK_ID - 1))
c_noise=${c_noise_values[$c_noise_index]}

# Run the Julia script with the current c_noise value
julia run_5.jl --T 2000 --c_noise $c_noise
