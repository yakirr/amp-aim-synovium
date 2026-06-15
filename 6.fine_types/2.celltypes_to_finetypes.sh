#!/bin/bash
function major_celltypes_to_finetypes()
{
    local ct=$1
    ct="${ct// /_}"
    local xen_path=$2
    local xen_path_escaped="${xen_path// /\\ }"
    chunk=$(echo "$xen_path" | grep -oP '_chunk\d+_' | grep -oP '\d+')
    echo $chunk
    if [ -z "$chunk" ]; then
    	echo "ERROR: could not extract chunk number from xen_path: $xen_path"
   	exit 1
    fi
    local cohort=$3
    local cca_path=$4
    local batch_vars=$5 
    
    filesize=$(stat -c%s "$xen_path")
    local mb=$((filesize / 1024 / 1024)) 
    if [ $mb -lt 200 ]; then
	mem="60G"
	time="00:45:00"
    elif [ $mb -lt 500 ]; then 
        mem="96G"
        time="01:00:00"
    elif [ $mb -lt 1000 ]; then 
	mem="200G"
	time="01:30:00"
    elif [ $mb -lt 2000 ]; then 
	mem="250G"
	time="02:00:00"
    else 
	mem="300G"
	time="03:00:00"
    fi
    
    if [ "$CONDA_DEFAULT_ENV" != "amp_harmony" ]; then
        echo "ERROR: wrong conda environment. Expected 'amp_harmony', got '$CONDA_DEFAULT_ENV'"
        exit 1
    fi
cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct}_${chunk}_majortypetofinetype.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct}_${chunk}_majortypetofinetype.err
#SBATCH -J ${cohort}_${ct}_${chunk}_majortypetofinetype
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
Rscript ./2.celltypes_to_finetypes.r "$ct" $xen_path_escaped "$cohort" "$chunk" "$cca_path" "$batch_vars" 
EOF
}
jq -r '.ct, .xen_path, .cohort, .cca_path, .batch_vars' "$1" \
| while read -r ct && read -r xen_path && read -r cohort && read -r cca_path && read -r batch_vars
do
    major_celltypes_to_finetypes "$ct" "$xen_path" "$cohort" "$cca_path" "$batch_vars"
done
