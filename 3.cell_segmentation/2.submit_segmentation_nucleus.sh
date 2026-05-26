#!/bin/bash
#usage: supply a list of sample ids via stdin (cat <file> | sh <this script>). this script
# will submit the baysor jobs corresponding to each sample id.

dir=/data/srlab/AMP_collab/lakshay-yakir/3.cell_segmentation
scripts=$dir/baysor

while IFS= read -r sid || [ -n "$sid" ]
do
    [ -z "$sid" ] && continue
    sample="$dir/out/$sid"
    echo "============== $sample ================"
    ${scripts}/baysor-subonesample-nucleus.sh "$sample/chunks/"
done
