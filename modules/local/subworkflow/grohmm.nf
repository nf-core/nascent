/*
 * TODO
 */

params.makeucscfile_options      = [:] // Collapses both strands, used as default value
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

include { GROHMM_MAKEUCSCFILE } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options  )
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
    // Generate UCSC files
    GROHMM_MAKEUCSCFILE ( bam )
    GROHMM_TRANSCRIPTCALLING_INDIVIDUAL ( bam )

    // Run Meta
    PICARD_MERGESAMFILES ( bam.collect() ).bam | GROHMM_TRANSCRIPTCALLING_META

    emit:
    transcripts      = GROHMM_TRANSCRIPTCALLING_INDIVIDUAL.out.transcripts
    meta_transcripts = GROHMM_TRANSCRIPTCALLING_META.out.transcripts
}
