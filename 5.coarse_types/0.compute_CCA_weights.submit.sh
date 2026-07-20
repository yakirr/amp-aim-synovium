#!/bin/bash
function compute_ct_cca()
{
    local ct=$1

    if [ "$CONDA_DEFAULT_ENV" != "amp_harmony" ]; then
   	echo "ERROR: wrong conda environment. Expected 'amp_harmony', got '$CONDA_DEFAULT_ENV'"
    	exit 1
    fi

cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/ 
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/compute_${ct}_cca.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/compute_${ct}_cca.err
#SBATCH -J compute_${ct}_cca
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH -c 1
#SBATCH -p normal,bigmem,long
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL
# ------------------------- End of Header ------------------------- #

#source activate /data/srlab/lsood/miniforge3/envs/amp_harmony 

Rscript ./0.compute_CCA_weights.r "$ct"

EOF
}

celltypes=("B" "Plasma" "T") 
for ct in ${celltypes[@]}; do 
	compute_ct_cca "$ct"
done

