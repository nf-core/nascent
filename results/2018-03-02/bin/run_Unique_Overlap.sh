#!/bin/bash

input_IMR90_data='../2018-02-06/bam'

# Sort GM eRNA
sortBed -i ../../data/2018-01-25/eRNA_GM_hg19.bed > eRNA_GM_hg19_sorted.bed                                                    

# IMR90 and overlap Analysis
mkdir IMR90_All

IMR90_Overlap='../2018-01-27/eRNA_IMR90_Overlap.bed'
output_IMR90_Overlap='./IMR90_All/'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting Overlap $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    coverageBed -d -sorted -a $i -b $IMR90_Overlap $IMR90_No_Overlap > ${i##*/};
    echo "Finished Overlap $i"
done ;

# Rename and Move everything to the right folder
mv *.bam $output_IMR90_Overlap
