#!/bin/bash

input_IMR90_data='../2018-02-06/bam'

# IMR90 Unique eRNAs Analysis
mkdir IMR90_Unique

IMR90_No_Overlap='../2018-01-29/eRNA_IMR90_hg19_No_Overlap.bed'
output_IMR90_Unique='./IMR90_Unique/'

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

# End IMR90 Unique Analysis 

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

# End IMR90 and overlap Analysis

# IMR90 and overlap Analysis
mkdir Overlap

Overlap='../2018-01-27/eRNA_IMR90_Overlap.bed'
output_Overlap='./IMR90_All/'

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

# End IMR90 and overlap Analysis

# GM and overlap Analysis
mkdir GM_Overlap

GM_Unique='../../data/2018-01-25/eRNA_GM_hg19.bed'
output_GM_Overlap='./GM_Overlap/'

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

# End IMR90 and overlap Analysis

# GM Unique Analysis
mkdir GM_Unique

output_GM_Unique='./GM_Unique/'

for i in `ls ${input_IMR90_data}/*.bam | cut -d "_" -f 1` ;
do
    echo "Starting GM $i"
    sort -k 1,1 -k2,2n $i
    # echo ${i##*/}
    coverageBed -d -sorted -a $i -b $GM_Unique > ${i##*/};
    echo "Finished GM $i"
done ;

# Rename and Move everything to the right folder
mv *.bam $output_GM_Unique

# End GM Unique Analysis
