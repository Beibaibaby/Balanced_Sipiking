#!/bin/bash

# Navigate to the directory containing the SLURM scripts
cd /home/dracoxu/Balanced_Sipiking

# Create a directory to store the job identifiers
mkdir -p id_store

# Generate a unique identifier based on the current date and time
JOB_ID=$(date '+%Y%m%d_%H%M%S')

# Write the unique identifier to a file in the id_store directory
echo $JOB_ID > id_store/unique_job_id

# Submit the SLURM job array and capture its job ID
array_job_id=$(sbatch pl.sbatch | cut -d ' ' -f 4)

# Prepare the SLURM script for processing the results
cat > process_results.sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=process_results
#SBATCH --account=dracoxu
#SBATCH --partition=tier1q
#SBATCH --time=01:00:00  # adjust time as needed
#SBATCH --mem=10gb       # adjust memory as needed
#SBATCH --dependency=afterok:$array_job_id

# Read the job identifier from the file
JOB_ID=\$(cat id_store/unique_job_id)

# Define the main directory where all the results are stored
main_dir="/gpfs/data/doiron-lab/draco/results_corr/exp_\${JOB_ID}"

# Run the Julia script to process and plot the results
julia plot_parallel.jl --main_dir "\$main_dir"
EOF

# Submit the processing job
sbatch process_results.sbatch

