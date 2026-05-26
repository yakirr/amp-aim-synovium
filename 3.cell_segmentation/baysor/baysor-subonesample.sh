chunks_folder=$1
xcoord=x
ycoord=y
zcoord=z
gene_col=gene
prior_seg=cell_ID

for tfile in $(find $chunks_folder -maxdepth 1 -iname "*.csv")
do
    echo $tfile

    # sample name = parent directory
    sample_name=$(basename "$(dirname "$(dirname "$tfile")")")
    echo "$sample_name"

    # chunk / FOV = filename without .csv
    sample_fov=$(basename "$tfile" .csv)
    echo "$sample_fov"
    
    # job name for Slurm
    job_name="B_${sample_name}_${sample_fov}"

    output_file="${chunks_folder}/../baysor_out/${sample_fov}"
    seg_file="${output_file}/segmentation.csv"
    if [[ -f "$seg_file" ]]; then
        echo -e "\e[32m$(basename "$seg_file") already exists\e[0m"
    else
        echo "$output_file"
        mkdir -p "$output_file"

        ntx=$(wc -l < "$tfile")
        echo "$ntx"
        ntx=$(awk -v val="$ntx" 'BEGIN {printf "%.2f", val/1000000}')
        if awk "BEGIN {exit !($ntx <= 0.5)}"; then
            mem_baysor=100G
            time_baysor="00-2:30"
            part_baysor=short
        elif awk "BEGIN {exit !($ntx <= 0.8)}"; then
            mem_baysor=200G
            time_baysor="00-5:00"
            part_baysor=normal
        elif awk "BEGIN {exit !($ntx <= 1.2)}"; then
            mem_baysor=200G
            time_baysor="00-8:00"
            part_baysor=normal
        elif awk "BEGIN {exit !($ntx <= 1.6)}"; then
            mem_baysor=200G
            time_baysor="00-11:30"
            part_baysor=normal
        else
            mem_baysor=250G
            time_baysor="00-12:00"
            part_baysor=normal
        fi

        cmd="sbatch --mem=$mem_baysor --time=$time_baysor --partition=$part_baysor \
            --job-name=$job_name \
            /data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation/baysor/baysor-run.sh \
            $output_file $tfile $xcoord $ycoord $zcoord $gene_col $prior_seg"
        echo -e "\e[33m$cmd\e[0m"
        eval "$cmd"
    fi
done
