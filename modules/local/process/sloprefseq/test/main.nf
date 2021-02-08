#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { BEDTOOLS_SLOPREFSEQ } from '../main.nf' addParams( options: [args: '-l 15 -r 15'] )

// Define input channels
// Run the workflow
workflow test_sloprefseq {
    def input = []
    input = [ [ id:'test'],
               file("${baseDir}/input/A.bed", checkIfExists: true),
                file("${baseDir}/input/genome.sizes", checkIfExists: true) ]

    BEDTOOLS_SLOPREFSEQ         ( input )


}
