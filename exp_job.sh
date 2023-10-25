#!/bin/bash -l

### SLURM Directives
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=jobname.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=jobname.%j.err  # Stderr (%j expands to jobId)

### Your experiment commands go here
julia run_inject_2.0.0.jl --T 50000 --f 0.0