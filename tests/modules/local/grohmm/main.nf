#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GROHMM_MAKEUCSCFILE } from '../../../../modules/local/process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_TRANSCRIPTCALLING } from '../../../../modules/local/process/grohmm/transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_PARAMETERTUNING } from '../../../../modules/local/process/grohmm/parametertuning/main.nf' addParams( options: [publish_dir:'test_grohmm'] )

workflow test_grohmm {
    def input = []
    input = [ [ id:'test'],
              [ file("${baseDir}/input/GroHMM_tutorial_input.bam", checkIfExists: true) ] ]
    GROHMM_MAKEUCSCFILE ( input )
}
