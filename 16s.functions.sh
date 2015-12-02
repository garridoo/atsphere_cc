#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

log() {

    echo $(date -u)": "$1 >> $logfile

}

cleanup() {
 
    find $working_dir -type f -print0 | xargs -0 rm -f

}

sampleSizes() {
    
    local seqs=$1
    local mapping_file=$2
    local output=$3

    rm -f $output

    for sample in $(cut -f 1 $mapping_file | tail -n +2)
    do 

        echo -n -e ""$sample"\t" >> $output
        echo $(grep -c ">"$sample"_" $seqs) >> $output

    done

    sort -n -k2,2 $output -o $output

}

