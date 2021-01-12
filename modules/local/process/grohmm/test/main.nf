#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GROHMM_MAKEUCSCFILE } from '../makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_TRANSCRIPTCALLING } from '../transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )

// Define input channels
// Run the workflow
workflow test_grohmm {
    def input = []
    input = [ [ id:'test'],
              [ file("${baseDir}/input/test.paired_end.name.sorted.bam", checkIfExists: true) ] ]
    GROHMM_MAKEUCSCFILE ( input )
}
