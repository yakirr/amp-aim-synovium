#!/bin/bash
function umap_full_ct_obj()
{
    local ct=$1
    ct="${ct// /_}"
    local lineage=$2 
    local xen_path=$3
    local xen_path_quoted="'${xen_path}'"
    local cohort=$4
    local cca_path=$5
    local batch_vars=$6
    filesize=$(stat -c%s "$xen_path")
    local mb=$((filesize / 1024 / 1024)) 
    if [ $mb -lt 2000 ]; then
	mem="250G"
	time="04:00:00"
    elif [ $mb -lt 5000 ]; then 
        mem="200G"
	time="06:00:00"
    elif [ $mb -lt 10000 ]; then 
	mem="300G"
	time="08:00:00"
    elif [ $mb -lt 20000 ]; then 
	mem="400G"
	time="10:00:00"
    else 
	mem="500G"
	time="12:00:00"
    fi
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
Rscript ./4.umap_combined_coarsetype_obj.r "$ct" "$lineage" $xen_path_quoted "$cohort" "$cca_path" "$batch_vars"
EOF
}
jq -r '.ct, .lineage, .xen_path, .cohort, .cca_path, .batch_vars' "$1" \
| while read -r ct && read -r lineage && read -r xen_path && read -r cohort && read -r cca_path && read -r batch_vars
do
    umap_full_ct_obj "$ct" "$lineage" "$xen_path" "$cohort" "$cca_path" "$batch_vars"
done
