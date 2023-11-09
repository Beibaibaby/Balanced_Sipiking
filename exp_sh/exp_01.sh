#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=jobname.%A_%a.out # %A is replaced by the job ID and %a by the array index
#SBATCH --error=jobname.%A_%a.err
#SBATCH --array=1-9%5 # Creates a job array with 9 tasks and a maximum of 5 running simultaneously

# Fixed parameters
d=0.2
f=0.95
T=3000
env=2
stim_start_time=400
stim_duration=5

# Define an array of stimstr values
stimstr_values=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

# Use SLURM_ARRAY_TASK_ID to select the stimstr value
stimstr=${stimstr_values[$SLURM_ARRAY_TASK_ID]}

# Run the Julia script with the selected stimstr value
julia run_inject_3.0.0.jl --d $d --f $f --T $T --env $env --stimstr $stimstr --stim_start_time $stim_start_time --stim_duration $stim_duration
