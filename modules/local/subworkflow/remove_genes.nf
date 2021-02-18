#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define input channels
// Run the workflow
include { BEDTOOLS_SORT } from '../../nf-core/software/bedtools/sort/main' addParams( options: [:] )
include { BEDTOOLS_SLOP } from '../../nf-core/software/bedtools/slop/main' addParams( options: [:] )
include { BEDTOOLS_INTERSECT } from '../../nf-core/software/bedtools/intersect/main' addParams( options: [:] )

workflow REMOVE_GENES {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    refseq // file: /path/to/refseq
    chromsizes

    main:
    BEDTOOLS_SORT ( refseq )
    BEDTOOLS_SLOP ( BEDTOOLS_SORT.out.bed, chromsizes )
    // BEDTOOLS_INTERSECT ( BEDTOOLS_SLOP.out.bed, reads )


    // emit:
    // transcrips = BEDTOOLS_INERSECT.out.bed
}
