#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HOMER_GROSEQ as HOMER_GROSEQ_BAM } from '../../../../../modules/nf-core/subworkflow/homer_groseq' addParams( options: [ publish_dir:'test_homer_groseq' ] )
include { HOMER_GROSEQ as HOMER_GROSEQ_BED } from '../../../../../modules/nf-core/subworkflow/homer_groseq' addParams( options: [ publish_dir:'test_homer_groseq' ] )

workflow test_homer_groseq_bam {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam', checkIfExists: true)]]
    def fasta = [ file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/genomics/sarscov2/fasta/test_genome.fasta', checkIfExists: true) ]

    HOMER_GROSEQ_BAM ( input, fasta )
}

workflow test_homer_groseq_bed {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/genomics/sarscov2/bed/test.bed', checkIfExists: true)]]
    def fasta = [ file('https://raw.githubusercontent.com/nf-core/modules/master/tests/data/genomics/sarscov2/fasta/test_genome.fasta', checkIfExists: true) ]

    HOMER_GROSEQ_BED ( input, fasta )
}
