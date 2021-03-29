/*
 * TODO
 */

params.makeucscfile_options      = [:] // Collapses both strands, used as default value
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

include { GROHMM_MAKEUCSCFILE      } from '../../modules/local/process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options  )
include { GROHMM_TRANSCRIPTCALLING } from '../../modules/local/process/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { PICARD_MERGESAMFILES     } from '../../modules/nf-core/software/picard/mergesamfiles/main'        addParams( options: [:]            )
include { GROHMM_PARAMETERTUNING   } from '../../modules/local/process/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


/*
 * Test with forward strand
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    PICARD_MERGESAMFILES ( bam.collect() )
    // Generate UCSC files
    GROHMM_MAKEUCSCFILE ( PICARD_MERGESAMFILES.out.bam )
    // Run Meta
    GROHMM_TRANSCRIPTCALLING ( PICARD_MERGESAMFILES.out.bam )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed         = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    wig         = GROHMM_MAKEUCSCFILE.out.wig
    plus_wig    = GROHMM_MAKEUCSCFILE.out.minuswig
    minus_wig   = GROHMM_MAKEUCSCFILE.out.pluswig
}
