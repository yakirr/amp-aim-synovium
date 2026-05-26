chunk_name=$1
xcoord=x
ycoord=y
zcoord=z
gene_col=gene
prior_seg=cell_ID

sample_name=${chunk_name/%"_"*}
sample_fov=${chunk_name/#${sample_name}_}

echo ${sample_name}
echo ${sample_fov}
chunks_folder="1.samples/${sample_name}/chunks"
output_file="${chunks_folder}/${sample_name}_${sample_fov}"
tfile="${chunks_folder}/${sample_name}_${sample_fov}.csv"
echo $output_file

mkdir -p $output_file
ntx=$(cat ${tfile} | wc -l)
echo ${ntx}
ntx=$(echo "scale=2; $ntx/1000000" | bc)
if (( $(echo "$ntx <= 0.5" | bc -l) )); then 
    mem_baysor=100G
    time_baysor="00-2:30"
    part_baysor=short
elif (( $(echo "$ntx <= 0.8" | bc -l) )); then
    mem_baysor=200G
    time_baysor="00-5:00"
    part_baysor=normal
elif (( $(echo "$ntx <= 1.2" | bc -l) )); then
    mem_baysor=200G
    time_baysor="00-8:00"
    part_baysor=normal
elif (( $(echo "$ntx <= 1.6" | bc -l) )); then
    mem_baysor=200G
    time_baysor="00-11:30"
    part_baysor=normal
else
    mem_baysor=250G
    time_baysor="00-12:00"
    part_baysor=normal
fi

cmd="sbatch --mem=$mem_baysor --time=$time_baysor --partition=$part_baysor \
    /data/srlab/AMP_collab/yakir/sh/baysor-run.sh \
    $output_file $tfile $xcoord $ycoord $zcoord $gene_col $prior_seg"
rm -f $output_file/out.txt
rm -f $output_file/err.txt
echo $cmd
eval "$cmd"
