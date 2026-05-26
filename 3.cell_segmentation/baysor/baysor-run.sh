#!/bin/bash

#SBATCH -c 1
#SBATCH -t 00-1:00
#SBATCH -p medium 
#SBATCH --mem 300G
#SBATCH --output=/data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation/slurm_logs/%x-%j.out
#SBATCH --error=/data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation/slurm_logs/%x-%j.err
output_file=$1
input_file=$2
xcoord=$3
ycoord=$4
zcoord=$5
gene_col=$6
prior_seg=$7

/data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation/baysor/baysor_package/bin/baysor run \
    -x $xcoord -y $ycoord -z $zcoord --gene-column $gene_col \
    --scale-std=100% --n-clusters=12 --prior-segmentation-confidence 0.7 \
    --polygon-format='FeatureCollection' --count-matrix-format='tsv' \
    -m 10 \
    -c '/data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation/baysor/xen_baysor_config.toml' \
    -o $output_file \
    $input_file :$prior_seg
     
