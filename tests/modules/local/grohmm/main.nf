#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GROHMM_MAKEUCSCFILE } from '../../../../modules/local/grohmm/makeucscfile/main.nf'
include { GROHMM_TRANSCRIPTCALLING } from '../../../../modules/local/grohmm/transcriptcalling/main.nf'
include { GROHMM_PARAMETERTUNING } from '../../../../modules/local/grohmm/parametertuning/main.nf'

workflow test_grohmm_makeucscfile {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['grohmm']['s0mR1'], checkIfExists: true),
            file(params.test_data['grohmm']['s40mR1'], checkIfExists: true) ]]
    GROHMM_MAKEUCSCFILE ( input )
}

workflow test_grohmm_transcriptcalling {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['grohmm']['s0mR1'], checkIfExists: true),
            file(params.test_data['grohmm']['s40mR1'], checkIfExists: true) ]]
    gtf = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)

    GROHMM_TRANSCRIPTCALLING (
        input,
        gtf,
        []
    )
}

workflow test_grohmm_parametertuning {
    def input = []
    input = [[ id: 'test' ],
            [ file(params.test_data['grohmm']['s0mR1'], checkIfExists: true),
            file(params.test_data['grohmm']['s40mR1'], checkIfExists: true) ]]
    GROHMM_PARAMETERTUNING ( input )
}
