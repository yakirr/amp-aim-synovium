function split_lineage_into_chunks()
{
    local lineage=$1
    local cohort=$2
    local max_mult=$3

    xen_path="/data/srlab/AMP_collab/lakshay-yakir/4.lineages/out_rds/${cohort}/lineages/${lineage}.rds"
    filesize=$(stat -c%s $xen_path)
    local mb=$((filesize / 1024 / 1024)) 

    if [ $mb -lt 200 ]; then
	mem="20G"
	time="00:30:00"
    elif [ $mb -lt 500 ]; then 
        mem="25G"
	time="00:45:00"
    elif [ $mb -lt 1000 ]; then 
	mem="30G"
	time="01:00:00"
    elif [ $mb -lt 2000 ]; then 
	mem="35G"
	time="01:30:00"
    else 
	mem="70G"
	time="02:00:00"
    fi
cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_splitchunks.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_splitchunks.err 
#SBATCH -J ${cohort}_${lineage}_splitchunks
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH -c 1
#SBATCH -p normal,bigmem,long
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL
# ------------------------- End of Header ------------------------- #

source activate /data/srlab/lsood/miniforge3/envs/amp_harmony 
echo "$lineage" 
echo "${mb} MB"
echo "$mem"
echo "$time"
Rscript ./1.split_lineages_into_chunks.r "$lineage" "$cohort" "$max_mult" 

EOF
}
jq -r '[.lineage, .cohort, .max_mult] | @tsv' "$1" \
| while IFS=$'\t' read -r lineage cohort max_mult
do 
    split_lineage_into_chunks $lineage $cohort $max_mult    
done
