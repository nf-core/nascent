#!/bin/bash

# A script to run bowtie2 on all files in a directory
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require Bowtie2 but it's not installed.  Aborting."; exit 1; }

# NOTE Correct the genome location on local computer
genomePath='/media/enhancer/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
fastqPath='/../../../data/2017-06-21/fastq'
# Make sure the Genome is there
if [ -f "$genomePath" ]
then
	echo "$genomePath found."
else
	echo "$genomePath not found."
fi

for i in `ls ${fastqPath}/*.fastq | cut -d "_" -f 1` ;
do
    echo "Start $i"
    bowtie2 --very-sensitive -t -p 4 -x $genomePath $i -S ${i##*/} ;
    echo "Finish $i"
done ;
