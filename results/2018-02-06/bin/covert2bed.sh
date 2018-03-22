#!/bin/bash

input='../../2017-06-27/sam/'

for file in ${input}/*.sam ;
    do sam2bed --keep-header < $file > ${file/%.sam/.bed};
done;

mkdir bed
mv *.bed bed/
