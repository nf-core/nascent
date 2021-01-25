#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GROHMM_MAKEUCSCFILE } from '../../../../modules/local/process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_TRANSCRIPTCALLING } from '../../../../modules/local/process/grohmm/transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_PARAMETERTUNING } from '../../../../modules/local/process/grohmm/parametertuning/main.nf' addParams( options: [publish_dir:'test_grohmm'] )

workflow test_grohmm_makeucscfile {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_MAKEUCSCFILE ( input )
}

workflow test_grohmm_transcriptcalling {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_TRANSCRIPTCALLING ( input )
}

workflow test_grohmm_parametertuning {
    def input = []
    input = [[ id: 'test' ],
             [ file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S0mR1.bam?raw=true', checkIfExists: true),
              file('https://github.com/Kraus-Lab/groHMM/blob/master/inst/extdata/S40mR1.bam?raw=true', checkIfExists: true) ]]
    GROHMM_PARAMETERTUNING ( input )
}
