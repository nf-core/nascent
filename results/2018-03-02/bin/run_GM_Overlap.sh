#!/bin/bash

input_IMR90_data='../../2018-02-06/bam'

# Sort GM eRNA
sortBed -i ../../../data/2018-01-25/eRNA_GM_hg19.bed > ../../../data/2018-01-25/eRNA_GM_hg19_sorted.bed                                                    

# GM and overlap Analysis
mkdir GM_Overlap

GM_Unique='../../../data/2018-01-25/eRNA_GM_hg19_sorted.bed'
output_GM_Overlap='../GM_Overlap/'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting GM and Overlap $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    coverageBed -d -sorted -a $i -b $GM_Unique $Overlap > ${i##*/};
    echo "Finished GM and Overlap $i"
done ;

# Rename and Move everything to the right folder
mv *.bam $output_GM_Overlap
