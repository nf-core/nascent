#!/bin/bash

# Make sure fastqc is installed
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require Fastqc but it's not installed.  Aborting."; exit 1; }

input_Fastq='../../data/2017-06-21/fastq'

for i in `ls ${input_Fastq}/*.fastq | cut -d "_" -f 1` ;
do
    echo "Starting $i"
    # echo ${i##*/}
    fastqc $i -o=fastqc/ -t=6
    echo "Finished $i"
done ;

mkdir fastqc && mv *.fastqc fastqc
