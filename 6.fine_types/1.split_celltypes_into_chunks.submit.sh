function split_ct_into_chunks()
{
    local ct="$1"
    local cohort="$2"
    local max_mult="$3"
    local ct_safe="${ct// /_}"   # for filenames/job names only

    xen_path="/data/srlab/AMP_collab/data/early_disease_synovium/xenium/combined/${cohort}/coarsetypes"
    filesize=$(stat -c%s "$xen_path")
    local mb=$((filesize / 1024 / 1024)) 
    if [ $mb -lt 200 ]; then
	mem="80G"
	time="00:30:00"
    elif [ $mb -lt 500 ]; then 
        mem="50G"
	time="00:45:00"
    elif [ $mb -lt 1000 ]; then 
	mem="60G"
	time="01:00:00"
    elif [ $mb -lt 2000 ]; then 
	mem="70G"
	time="01:30:00"
    else 
	mem="100G"
	time="02:00:00"
    fi
cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/6.fine_types
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct_safe}_splitchunks.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/6.fine_types/slurm_logs/${cohort}_${ct_safe}_splitchunks.err 
#SBATCH -J ${cohort}_${ct_safe}_splitchunks
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
Rscript ./1.split_celltypes_into_chunks.r "$ct" "$cohort" "$max_mult" 
EOF
}
jq -r '[.ct, .cohort, .max_mult] | @tsv' "$1" \
| while IFS=$'\t' read -r ct cohort max_mult
do 
    split_ct_into_chunks "$ct" "$cohort" "$max_mult"    
done
