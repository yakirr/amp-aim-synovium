#!/bin/bash
function compute_top50genes_overlap()
{
    local ct=$1
    local cohort=$2
    
    obj_dir="/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/out_rds/${cohort}/${ct}"
    filesize=$(du -sb "$obj_dir" | cut -f1)
    local mb=$((filesize / 1024 / 1024)) 

    if [ $mb -lt 10000 ]; then
	mem="40G"
	time="00:30:00"
    elif [ $mb -lt 15000 ]; then 
        mem="60G"
	time="01:00:00"
    elif [ $mb -lt 20000 ]; then 
	mem="80G"
	time="01:00:00"
    elif [ $mb -lt 30000 ]; then 
	mem="100G"
	time="01:00:00"
    else 
	mem="120G"
	time="02:00:00"
    fi

    if [ "$CONDA_DEFAULT_ENV" != "amp_harmony" ]; then
   	echo "ERROR: wrong conda environment. Expected 'amp_harmony', got '$CONDA_DEFAULT_ENV'"
    	exit 1
    fi
    
    ct="${ct// /_}"

cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/ 
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct}_compute_top50genes.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct}_compute_top50genes.err 
#SBATCH -J ${cohort}_${ct}_compute_top50genes
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH -c 1
#SBATCH -p normal,bigmem,long
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL
# ------------------------- End of Header ------------------------- #

#source activate /data/srlab/lsood/miniforge3/envs/amp_harmony 
echo "$ct" 
echo "${mb} MB"
echo "$mem"
echo "$time"

Rscript ./4.compute_top50genes_overlap.r "$ct" "EDP1-EDP2-ARB" 

EOF
}

#cts=("B" "Dendritic cell" "Lining" "Macrophage" "NK" "Plasma" "Sublining" "T" "Vascular endothelial")
cts=("Vascular endothelial" "Lining" "Macrophage" "Sublining")

for ct in "${cts[@]}"; do  
	compute_top50genes_overlap "${ct}" "EDP1-EDP2-ARB"
done 

