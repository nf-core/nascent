#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOPREFSEQ } from '../../../../modules/local/process/sloprefseq/main.nf' addParams( options: [args: '-l 15 -r 15'] )

// Define input channels
// Run the workflow
workflow test_sloprefseq {
    def input = []
    input = [ [ id:'test'],
                file("https://github.com/nf-core/modules/blob/master/tests/data/bed/A.bed", checkIfExists: true),
                file("https://github.com/nf-core/modules/blob/master/tests/data/bed/genome.sizes", checkIfExists: true) ]

    BEDTOOLS_SLOPREFSEQ         ( input )


}
