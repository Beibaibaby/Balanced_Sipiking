#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/out/demo_job.%j.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/out/demo_job.%j.err

# Run the Julia script
julia run_7_t150.jl --T 2000 --dir_name_in “/gpfs/data/doiron-lab/draco/results_corr/” --ie_sign true --ee_sign true --add_noise false --stim_duration 80 --stimstr 0.0
#5.4