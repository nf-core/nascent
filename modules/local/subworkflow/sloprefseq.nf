#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_SLOPREFSEQ } from '../process/sloprefseq/main.nf' addParams( options: [args: '-l 15 -r 30'] )

// Define input channels
// Run the workflow
workflow test_sloprefseq {
    def input = []
    input = [ [ id:'test'],
               //file("${baseDir}/input/A.bed", checkIfExists: true),
                //file("${baseDir}/input/genome.sizes", checkIfExists: true) ]
                file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/A.bed", checkIfExists: true),
                file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/genome.sizes", checkIfExists: true) ]

    BEDTOOLS_SLOPREFSEQ         ( input )


}

