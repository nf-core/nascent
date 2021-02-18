#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {BEDTOOLS_SORT } from '../../../../../modules/nf-core/software/bedtools/sort/main' addParams( options: [publish_dir:'test_removegenes'] )
include {BEDTOOLS_SLOP } from '../../../../../modules/nf-core/software/bedtools/slop/main' addParams( options: [publish_dir:'test_removegenes', args: '-l 15 -r 30'] )
include {BEDTOOLS_INTERSECT } from '../../../../../modules/nf-core/software/bedtools/intersect/main' addParams( options: [publish_dir:'test_removegenes'] )

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

workflow REMOVE_GENES {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    refseq // file: /path/to/refseq
    chromsizes

    main:
    BEDTOOLS_SORT ( refseq )
    BEDTOOLS_SLOP ( BEDTOOLS_SORT.out.bed, chromsizes )
}

