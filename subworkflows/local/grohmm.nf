/*
 * TODO
 */

params.makeucscfile_options      = [:] // Collapses both strands, used as default value
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

include { GROHMM_MAKEUCSCFILE      } from '../../modules/local/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options  )
include { GROHMM_TRANSCRIPTCALLING } from '../../modules/local/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { PICARD_MERGESAMFILES     } from '../../modules/nf-core/software/picard/mergesamfiles/main'        addParams( options: [:]            )
include { GROHMM_PARAMETERTUNING   } from '../../modules/local/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


/*
 * Test with forward strand
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    // FIXME
    // ch_meta = Channel.value([[id: 'meta', single_end:false ]])
    // meta_input = ch_meta.concat(  bam.collect{it[1]}.flatten() )
    // PICARD_MERGESAMFILES ( meta_input )

    // Generate UCSC files
    GROHMM_MAKEUCSCFILE ( bam )
    // Run Meta
    GROHMM_TRANSCRIPTCALLING ( bam )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed         = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    wig         = GROHMM_MAKEUCSCFILE.out.wig
    plus_wig    = GROHMM_MAKEUCSCFILE.out.minuswig
    minus_wig   = GROHMM_MAKEUCSCFILE.out.pluswig
}
