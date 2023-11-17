#!/bin/bash
#SBATCH --job-name=draco_test
#SBATCH --account=dracoxu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:30:30
#SBATCH --cpus-per-task=2
#SBATCH --mem=30gb
#SBATCH --output=/gpfs/data/doiron-lab/draco/trash/demo_job.%j.out
#SBATCH --error=/gpfs/data/doiron-lab/draco/trash/demo_job.%j.err

echo hello > /gpfs/data/doiron-lab/foo