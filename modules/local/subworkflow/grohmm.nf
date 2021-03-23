/*
 * TODO
 */

params.makeucscfile_options      = [:] // Collapses both strands, used as default value
params.transcriptcalling_options = [:]
params.parametertuning_options   = [:]

include { GROHMM_MAKEUCSCFILE      } from '../process/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options  )
include { GROHMM_TRANSCRIPTCALLING } from '../process/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { PICARD_MERGESAMFILES     } from '../../nf-core/software/picard/mergesamfiles/main'        addParams( options: [:]            )
include { GROHMM_PARAMETERTUNING   } from '../process/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


/*
 * Test with forward strand
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    // Generate UCSC files
    PICARD_MERGESAMFILES ( bam.collect() ).bam | GROHMM_MAKEUCSCFILE
    // Run Meta
    PICARD_MERGESAMFILES ( bam.collect() ).bam | GROHMM_TRANSCRIPTCALLING_META

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    wig         = GROHMM_MAKEUCSCFILE.out.wig
    plus_wig    = GROHMM_MAKEUCSCFILE.out.minuswig
    minus_wig   = GROHMM_MAKEUCSCFILE.out.pluswig
}
