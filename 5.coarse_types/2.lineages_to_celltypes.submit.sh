#!/bin/bash
function lineage_to_major_celltypes()
{
    local lineage=$1
    local xen_path=$2
    chunk=$(echo "$xen_path" | grep -oP 'chunk\d+' | grep -oP '\d+')
    if [ -z "$chunk" ]; then
    	echo "ERROR: could not extract chunk number from xen_path: $xen_path"
   	 exit 1
    fi
    local cohort=$3
    local batch_vars=$4 

    filesize=$(stat -c%s $xen_path)
    local mb=$((filesize / 1024 / 1024)) 

    if [ $mb -lt 200 ]; then
	mem="60G"
	time="01:30:00"
    elif [ $mb -lt 500 ]; then 
        mem="96G"
	time="03:00:00"
    elif [ $mb -lt 1000 ]; then 
	mem="200G"
	time="06:00:00"
    elif [ $mb -lt 2000 ]; then 
	mem="250G"
	time="09:00:00"
    else 
	mem="300G"
	time="12:00:00"
    fi
cat << EOF | sbatch
#!/bin/bash
#SBATCH -D /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/ 
#SBATCH -o /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_${chunk}_lineagetomajortype.out
#SBATCH -e /data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/slurm_logs/${cohort}_${lineage}_${chunk}_lineagetomajortype.err 
#SBATCH -J ${cohort}_${lineage}_${chunk}_lineagetomajortype
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
Rscript ./2.lineages_to_celltypes.r "$lineage" "$xen_path" "$cohort" "${chunk}" "$batch_vars" 

EOF
}
jq -r '[.lineage, .xen_path, .cohort, .batch_vars] | @tsv' "$1" \
| while IFS=$'\t' read -r lineage xen_path cohort batch_vars 
do
    #if [ ! -f "${outpath}" ]; then  
        lineage_to_major_celltypes "$lineage" "$xen_path" "$cohort" "$batch_vars" 
    #else 
    #	echo "Skipping ${ct} pipeline ${out} thresh ${gene_filter_thresh} iter ${iter} - output already exists" 
    #fi
done
