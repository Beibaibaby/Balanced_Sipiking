#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/demo_job.%j.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/demo_job.%j.err

# Run the Julia script
julia run_6.jl --T 200 --c_noise 0.5 --scale_noise 2.5 --dir_name_in "/home/dracoxu/Balanced_Sipiking"