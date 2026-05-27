#!/usr/bin/env bash
#SBATCH --job-name=cells_to_level3
#SBATCH --output=slurm_logs/cells_to_level3_%A_%a.out
#SBATCH --error=slurm_logs/cells_to_level3_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --partition=long,bigmem,normal
#SBATCH --cpus-per-task=1
#SBATCH --array=1-102

set -euo pipefail

mkdir -p slurm_logs

echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

echo "Loading R env"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate r-seurat-h5

echo "running script"
Rscript 4.cells_to_level3.R "${SLURM_ARRAY_TASK_ID}"

echo "Job finished: $(date)"