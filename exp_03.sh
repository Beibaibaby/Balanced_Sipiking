#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=jobname.%j.out
#SBATCH --error=jobname.%j.err

# Fixed parameters
d=0.01
f=0.98
T=3000
env=2
stim_start_time=400
stim_duration=5

# Loop over stimstr values
for stimstr in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
    julia run_inject_3.0.0.jl --d $d --f $f --T $T --env $env --stimstr $stimstr --stim_start_time $stim_start_time --stim_duration $stim_duration
done