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
#SBATCH --array=1-6%3 # Creates a job array with 6 tasks and a maximum of 3 running simultaneously

# Fixed parameters
f=0.98
T=3000
env=2
stimstr=0.8
stim_start_time=400
stim_duration=5

# Array of d values
d_values=(0 0.01 0.05 0.10 0.2 0.3 0.4)

# Use SLURM_ARRAY_TASK_ID to select the d value
d=${d_values[$SLURM_ARRAY_TASK_ID]}

# Run the Julia script with the selected d value
julia run_inject_3.0.0.jl --d $d --f $f --T $T --env $env --stimstr $stimstr --stim_start_time $stim_start_time --stim_duration $stim_duration
