#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { REMOVE_GENES } from '../../../../../modules/local/subworkflow/remove_genes.nf' addParams( options: [publish_dir:'test_removegenes'] )

// Define input channels
// Run the workflow
workflow test_remove_genes {
    def input = []
    input = [ [ id:'test'],
               file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/A.bed", checkIfExists: true),
                file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/B.bed", checkIfExists: true),
                file("https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/genome.sizes", checkIfExists: true) ]

    REMOVE_GENES         ( input )
}

