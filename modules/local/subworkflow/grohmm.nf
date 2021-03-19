/*
 * TODO
 */

params.makeucscfile_options      = [:]
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

// FIXME include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_FORWARD } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options, args:'--strand=+' )
// FIXME include { GROHMM_MAKEUCSCFILE as GROHMM_MAKEUCSCFILE_REVERSE } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile, args:'--strand=-' )
include {
    GROHMM_TRANSCRIPTCALLING as GROHMM_TRANSCRIPTCALLING_INDIVIDUAL
    GROHMM_TRANSCRIPTCALLING as GROHMM_TRANSCRIPTCALLING_META } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { PICARD_MERGESAMFILES                                } from '../../nf-core/software/picard/mergesamfiles/main'        addParams( options: [:]            )
include { GROHMM_PARAMETERTUNING                              } from '../process/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


/*
 * Test with forward strand
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    GROHMM_TRANSCRIPTCALLING_INDIVIDUAL ( bam )

    // Run Meta
    PICARD_MERGESAMFILES ( bam.collect() ).bam | GROHMM_TRANSCRIPTCALLING_META

    emit:
    transcripts      = GROHMM_TRANSCRIPTCALLING_INDIVIDUAL.out.transcripts
    meta_transcripts = GROHMM_TRANSCRIPTCALLING_META.out.transcripts
}
