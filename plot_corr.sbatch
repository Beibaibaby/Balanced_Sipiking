#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=22:30:30
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/a1.%A_%a.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/a1.%A_%a.err
#SBATCH --array=1-10%10

# Use SLURM_ARRAY_JOB_ID as the identifier
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}

export QT_QPA_PLATFORM='offscreen'

# Generate a directory name using the array job ID and task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_corr/exp_${ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
julia run_6.jl --T 3000 --dir_name_in "$dir_name_in"

# Create a file indicating the completion of this particular task

main_dir="/gpfs/data/doiron-lab/draco/results_corr/exp_${ARRAY_JOB_ID}"

julia plot_parallel.jl --main_dir "$main_dir"