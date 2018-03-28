#!/bin/bash


command -v sam2bed >/dev/null 2>&1 || { echo >&2 "I require bedops but it's not installed.  Aborting."; exit 1; }

input='../2017-06-27/sam'

for file in ${input}/*.sam ;
do
    echo "start $file"
    sam2bed --keep-header -i sam $file > ${file/%.sam/.bed};
    echo "finish $file"
done;

mkdir bed
mv ../2017-06-27/sam/*.bed bed/
