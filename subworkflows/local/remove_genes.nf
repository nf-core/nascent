#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.sort_options    = [:]
params.slop_options    = [:]
params.intersect_options = [:]

include { BEDTOOLS_SORT } from '../../modules/nf-core/modules/bedtools/sort/main' addParams( options: params.sort_options )
include { BEDTOOLS_SLOP } from '../../modules/nf-core/modules/bedtools/slop/main' addParams( options: params.slop_options )
include { BEDTOOLS_INTERSECT } from '../../modules/nf-core/modules/bedtools/intersect/main' addParams( options: params.intersect_options )

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
