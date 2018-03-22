#!/usr/bin/env bash

# Coverts sam to bam
input='../2017-06-27/sam'

for file in ${input}/*.sam ;
do
    echo "start $file"
    touch ${file##*/}
    samtools view -bS $file > ${file##*/};
    echo "finish $file"
done;

mkdir bam
mv *.sam *.bam
mv *.bam bam/
