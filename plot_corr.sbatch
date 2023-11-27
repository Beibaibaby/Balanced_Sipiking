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
#SBATCH --array=1-400%100

# Use SLURM_ARRAY_JOB_ID as the identifier
ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}

export QT_QPA_PLATFORM='offscreen'

# Generate a directory name using the array job ID and task ID
dir_name_in="/gpfs/data/doiron-lab/draco/results_corr/exp_${ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"

event_thre = 1.8

# Create the directory along with any necessary parent directories
mkdir -p "$dir_name_in"

# Run the Julia script with the current values and pass the generated directory name
julia run_6.jl --T 10000 --dir_name_in "$dir_name_in" --c_noise 0.1 --sigma_noise 8.0 --event_thre $event_thre

# Create a file indicating the completion of this particular task

main_dir="/gpfs/data/doiron-lab/draco/results_corr/exp_${ARRAY_JOB_ID}"

julia plot_parallel.jl --main_dir "$main_dir" --event_thre $event_thre