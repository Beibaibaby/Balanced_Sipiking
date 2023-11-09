#!/bin/bash -l

### SLURM Directives
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:50:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --array=1-7%7
#SBATCH --output=draco_test_%A_%a.out # Stdout (%A expands to job ID, %a to array index)
#SBATCH --error=draco_test_%A_%a.err  # Stderr (%A expands to job ID, %a to array index)

# Define an array of stimstr values
stimstr_values=(0.7 0.75 0.80)

# Get the stimstr value for the current array task
s=${stimstr_values[$SLURM_ARRAY_TASK_ID]}

# Load Julia module if required
# module load julia/1.x.x

# Execute the Julia script with the specified parameters and the current stimstr value
julia run_inject_3.0.0.jl --d 0.01 --f 0.98 --T 1000 --env 2 --stimstr $s --stim_start_time 400 --stim_duration 5
