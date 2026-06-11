#!/bin/bash
function umap_full_ct_obj()
{
    local ct=$1
    local lineage=$2 
    local cohort=$3
    local cca_path=$4
    local batch_vars=$5
    
    filesize=$(find "/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds/${cohort}/${lineage}/${ct}" \
    -name "*labeltransfer.rds" -print0 \
    | xargs -0 stat -c%s \
    | awk '{sum += $1} END {print sum+0}')
    local mb=$((filesize / 1024 / 1024)) 
    if [ $mb -lt 200 ]; then
	mem="30G"
	time="00:30:00"
    elif [ $mb -lt 500 ]; then 
        mem="50G"
	time="01:00:00"
    elif [ $mb -lt 1000 ]; then 
	mem="100G"
	time="01:30:00"
    elif [ $mb -lt 2000 ]; then 
	mem="200G"
	time="04:00:00"
    else 
	mem="250G"
	time="06:00:00"
    fi
    ct="${ct// /_}"
    if [ "$CONDA_DEFAULT_ENV" != "amp_harmony" ]; then
   	echo "ERROR: wrong conda environment. Expected 'amp_harmony', got '$CONDA_DEFAULT_ENV'"
    	exit 1
    fi
cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/ 
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${ct}_umap_combinedcelltype.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${ct}_umap_combinedcelltype.err 
#SBATCH -J ${cohort}_${ct}_umap_combinedcelltype
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
Rscript ./3.combine_coarsetype_chunks_and_umap.r "$ct" "$lineage" "$cohort" "$cca_path" "$batch_vars"
EOF
}
jq -r '.ct, .lineage, .cohort, .cca_path, .batch_vars' "$1" \
| while read -r ct && read -r lineage && read -r cohort && read -r cca_path && read -r batch_vars
do
    umap_full_ct_obj "$ct" "$lineage" "$cohort" "$cca_path" "$batch_vars"
done
