/*
 * TODO
 */

include { GROHMM_MAKEUCSCFILE      } from '../../modules/local/grohmm/makeucscfile/main.nf'
include { GROHMM_TRANSCRIPTCALLING } from '../../modules/local/grohmm/transcriptcalling/main.nf'
include { GROHMM_PARAMETERTUNING   } from '../../modules/local/grohmm/parametertuning/main.nf'

/*
 * Note meta refers to all merged files
 */
workflow GROHMM {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    gtf

    main:
    // Generate UCSC files
    GROHMM_MAKEUCSCFILE ( bam )

    ch_tuning = Channel.empty()

    if (params.with_tuning) {
        GROHMM_PARAMETERTUNING (bam, file("${launchDir}/modules/local/grohmm/test/tune.csv") )
        ch_tuning = GROHMM_PARAMETERTUNING.out.tuning
    } else {
        ch_tuning = []
    }

    GROHMM_TRANSCRIPTCALLING (
        bam,
        gtf,
        ch_tuning
    )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed         = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    td_plot     = GROHMM_TRANSCRIPTCALLING.out.td_plot
    wig         = GROHMM_MAKEUCSCFILE.out.wig
    plus_wig    = GROHMM_MAKEUCSCFILE.out.minuswig
    minus_wig   = GROHMM_MAKEUCSCFILE.out.pluswig
}
