/*
 * TODO
 */

params.makeucscfile_options         = [:] // Collapses both strands, used as default value
params.transcriptcalling_options    = [:]
params.parametertuning_options      = [:]

include { GROHMM_MAKEUCSCFILE      } from '../../modules/local/grohmm/makeucscfile/main.nf'      addParams( options: params.makeucscfile_options  )
include { GROHMM_TRANSCRIPTCALLING } from '../../modules/local/grohmm/transcriptcalling/main.nf' addParams( options: params.transcriptcalling_options )
include { GROHMM_PARAMETERTUNING   } from '../../modules/local/grohmm/parametertuning/main.nf'   addParams( options: params.parametertuning_options )


/*
 * Note meta refers to all merged files
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    // Generate UCSC files
    GROHMM_MAKEUCSCFILE ( bam )
    // Run Meta
   // GROHMM_PARAMETERTUNING (bam, file("${launchDir}/modules/local/grohmm/test/tune.csv") )
    GROHMM_TRANSCRIPTCALLING ( bam, [] )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed         = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    wig         = GROHMM_MAKEUCSCFILE.out.wig
    plus_wig    = GROHMM_MAKEUCSCFILE.out.minuswig
    minus_wig   = GROHMM_MAKEUCSCFILE.out.pluswig
}

