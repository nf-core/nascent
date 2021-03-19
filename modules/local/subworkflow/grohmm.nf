/*
 * TODO
 */

params.makeucscfile_options      = [:]
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

// FIXME include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_FORWARD } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options, args:'--strand=+' )
// FIXME include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_REVERSE } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile, args:'--strand=-' )
include { GROHMM_TRANSCRIPTCALLING                           } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { GROHMM_PARAMETERTUNING                             } from '../process/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


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
