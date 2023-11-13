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
julia run_5.jl --T 1500 --c_noise 0.0 --ie_sign false --ee_sign false