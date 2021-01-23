#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GROHMM_MAKEUCSCFILE } from '../process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm', norm:2] )
include { GROHMM_TRANSCRIPTCALLING } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_PARAMETERTUNING } from '../process/grohmm/parametertuning/main.nf' addParams( options: [publish_dir:'test_grohmm'] )

// Define input channels
// Run the workflow
workflow GROHMM {
    def input = []
    input = [ file("${baseDir}/modules/local/process/grohmm/test/input/S0mR1.bam", checkIfExists: true) ]
    GROHMM_PARAMETERTUNING ( input )
}
