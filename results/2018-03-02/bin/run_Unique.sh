#!/bin/bash

input_IMR90_data='../2018-02-06/bam'

# IMR90 Unique eRNAs Analysis
mkdir IMR90_Unique

IMR90_No_Overlap='../../2018-01-29/eRNA_IMR90_hg19_No_Overlap.bed'
output_IMR90_Unique='../IMR90_Unique/'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting Unique $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    coverageBed -d -sorted -a $i -b $IMR90_No_Overlap > ${i##*/};
    echo "Finished Unique $i"
done ;

# Rename and Move everything to the right folder
# for f in *.bam; do printf '%s\n' "${f%.bam}_Unique.bam"; done
mv *.bam $output_IMR90_Unique
