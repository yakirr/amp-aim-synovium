#!/bin/bash
function compute_top50genes_overlap()
{
    local lineage=$1
    local cohort=$2

    obj_dir="/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds/${cohort}/${lineage}"
    filesize=$(du -sb $obj_dir | cut -f1)
    local mb=$((filesize / 1024 / 1024)) 

    if [ $mb -lt 2000 ]; then
	mem="30G"
	time="00:30:00"
    elif [ $mb -lt 5000 ]; then 
        mem="40G"
	time="01:00:00"
    elif [ $mb -lt 10000 ]; then 
	mem="50G"
	time="01:00:00"
    elif [ $mb -lt 20000 ]; then 
	mem="75G"
	time="01:00:00"
    else 
	mem="100G"
	time="02:00:00"
    fi

    if [ "$CONDA_DEFAULT_ENV" != "amp_harmony" ]; then
   	echo "ERROR: wrong conda environment. Expected 'amp_harmony', got '$CONDA_DEFAULT_ENV'"
    	exit 1
    fi

cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/ 
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_compute_top50genes.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_compute_top50genes.err 
#SBATCH -J ${cohort}_${lineage}_compute_top50genes
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH -c 1
#SBATCH -p normal,bigmem,long
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL
# ------------------------- End of Header ------------------------- #

#source activate /data/srlab/lsood/miniforge3/envs/amp_harmony 
echo "$lineage" 
echo "${mb} MB"
echo "$mem"
echo "$time"

Rscript ./4.compute_top50genes_overlap.r "$lineage" "EDP1-EDP2-ARB" 

EOF
}

#lineages=("B_plasma" "Myeloid" "Stromal" "T_NK" "Endothelial")
lineages=("T_NK")
for lineage in ${lineages[@]}; do 
	compute_top50genes_overlap "$lineage" "EDP1-EDP2-ARB"
done

