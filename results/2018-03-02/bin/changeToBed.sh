#!/bin/bash

input_IMR90_data='../../2018-02-06/bam'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    bamToBed -i $i > ${i##*/};
    echo "Finished $i"
done ;

mv *.bam *.bed
mkdir ../input_IMR90/
mv ./*.bed ../input_IMR90/
