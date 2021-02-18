#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.remove_genes_options = [ publish_dir: 'test_removegenes', slop_options: '-l 15 -r 30' ]

include { REMOVE_GENES } from '../../../../../modules/local/subworkflow/remove_genes.nf' addParams( options: params.remove_genes_options )

// Define input channels
// Run the workflow
workflow test_remove_genes {
    def reads = []
    def refseq = []
    def chromsizes  = []
    reads = [ [ id:'test'],
               file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/A.bed', checkIfExists: true) ]
    refseq = [ [ id:'test1'],
               file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/B.bed', checkIfExists: true) ]
    chromsizes = [ file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/bed/genome.sizes', checkIfExists: true) ]

    REMOVE_GENES         ( reads,refseq,chromsizes )
}
