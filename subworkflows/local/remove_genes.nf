#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SORT } from '../../modules/nf-core/modules/bedtools/sort/main'
include { BEDTOOLS_SLOP } from '../../modules/nf-core/modules/bedtools/slop/main'
include { BEDTOOLS_INTERSECT } from '../../modules/nf-core/modules/bedtools/intersect/main'

workflow REMOVE_GENES {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    refseq // file: /path/to/refseq
    chromsizes

    main:
    BEDTOOLS_SORT ( refseq )
    BEDTOOLS_SLOP ( BEDTOOLS_SORT.out.bed, chromsizes )
    BEDTOOLS_INTERSECT ( BEDTOOLS_SLOP.out.bed, reads )


    emit:
    transcripts = BEDTOOLS_INTERSECT.out.bed
}
