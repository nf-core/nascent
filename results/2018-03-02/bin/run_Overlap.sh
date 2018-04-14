#!/bin/bash

input_IMR90_data='../../2018-02-06/bam'

# overlap Coverage
mkdir Overlap

Overlap='../../2018-01-27/eRNA_IMR90_Overlap.bed'
output_Overlap='../Overlap/'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting Overlap $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    coverageBed -d -sorted -a $i -b $Overlap > ${i##*/};
    echo "Finished Overlap $i"
done ;

# Rename and Move everything to the right folder
mv *.bam $output_Overlap
