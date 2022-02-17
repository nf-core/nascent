#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_GROSEQ as HOMER_GROSEQ_BAM } from '../../../../subworkflows/nf-core/homer_groseq'
include { HOMER_GROSEQ as HOMER_GROSEQ_BED } from '../../../../subworkflows/nf-core/homer_groseq'

workflow test_homer_groseq_bam {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['nf-core']['test_paired_end_sorted_bam'], checkIfExists: true)]]
    def fasta = [ file(params.test_data['nf-core']['genome_fasta'], checkIfExists: true) ]

    HOMER_GROSEQ_BAM ( input, fasta )
}

workflow test_homer_groseq_bed {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['nf-core']['bed']['test_bed'], checkIfExists: true)]]
    def fasta = [ file(params.test_data['nf-core']['genome_fasta'], checkIfExists: true) ]

    HOMER_GROSEQ_BED ( input, fasta )
}
