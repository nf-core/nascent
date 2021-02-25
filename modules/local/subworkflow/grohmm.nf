#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_FORWARD} from '../process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm', args:'--strand=+'] )
include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_REVERSE } from '../process/grohmm/makeucscfile/main.nf' addParams( options: [publish_dir:'test_grohmm', args:'--strand=-'] )
include { GROHMM_TRANSCRIPTCALLING } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: [publish_dir:'test_grohmm'] )
include { GROHMM_PARAMETERTUNING } from '../process/grohmm/parametertuning/main.nf' addParams( options: [publish_dir:'test_grohmm'] )


/*
 * Test with forward strand
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    GROHMM_TRANSCRIPTCALLING ( bam )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
}
