#!/bin/bash
#SBATCH --job-name=combine_cells
#SBATCH --output=slurm_logs/combine_cells_%j.out
#SBATCH --error=slurm_logs/combine_cells_%j.err
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

# Create log directory if it doesn't exist
mkdir -p slurm_logs

echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

echo "Loading R env"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate r2026

echo "running script"
Rscript 3a.combine_cells.R

echo "Job finished: $(date)"