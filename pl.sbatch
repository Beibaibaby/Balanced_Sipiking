#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:30:30
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/aa.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/aa.%A_%a.err
#SBATCH --array=1-10%30

# Read the common job identifier from the file
JOB_ID=$(cat /home/dracoxu/Balanced_Sipiking/unique_job_id)

# Generate a directory name using the common job identifier and array task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_corr/exp_${JOB_ID}/${SLURM_ARRAY_TASK_ID}"

# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
julia run_5_parallel.jl --T 200 --dir_name_in "$dir_name_in"

# Check if the current task is the last one in the array
if [ "$SLURM_ARRAY_TASK_ID" -eq 10 ]; then
    # Define the main directory where all the results are stored
    main_dir="/gpfs/data/doiron-lab/draco_results/corr/exp_${JOB_ID}"

    # Run the Julia script to process and plot the results
    julia plot_parallel.jl --main_dir "$main_dir"
fi
