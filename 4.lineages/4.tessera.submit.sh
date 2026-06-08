#!/bin/bash
#SBATCH --job-name=tessera
#SBATCH --output=slurm_logs/4.tessera_%j.out
#SBATCH --error=slurm_logs/4.tessera_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal,bigmem,long

set -euo pipefail

# Create log directory if it doesn't exist
mkdir -p slurm_logs

echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

echo "Loading R env"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate r2026

echo "running script"
Rscript 4.tessera.r

echo "Job finished: $(date)"