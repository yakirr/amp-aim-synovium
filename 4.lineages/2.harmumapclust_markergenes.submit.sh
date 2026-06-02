#!/bin/bash
#SBATCH --job-name=harmumapclust
#SBATCH --output=slurm_logs/2.harmumapclust_markergenes_%j.out
#SBATCH --error=slurm_logs/2.harmumapclust_markergenes_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=999G
#SBATCH --cpus-per-task=16
#SBATCH --partition=bigmem

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
Rscript 2.harmumapclust_markergenes.r

echo "Job finished: $(date)"