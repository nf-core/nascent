/*
 * TODO
 */

include { GROHMM_TRANSCRIPTCALLING } from '../../modules/local/grohmm/transcriptcalling/main.nf'
include { GROHMM_PARAMETERTUNING   } from '../../modules/local/grohmm/parametertuning/main.nf'

/*
 * Note meta refers to all merged files
 */
workflow GROHMM {
    take:
    bams // channel: [ val(meta), [ bams ] ]
    gtf

    main:
    ch_tuning = Channel.empty()

    if (params.with_tuning) {
        GROHMM_PARAMETERTUNING (bam, file("${launchDir}/modules/local/grohmm/test/tune.csv") )
        ch_tuning = GROHMM_PARAMETERTUNING.out.tuning
    } else {
        ch_tuning = []
    }

    GROHMM_TRANSCRIPTCALLING (
        bams,
        gtf,
        ch_tuning
    )

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed         = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    td_plot     = GROHMM_TRANSCRIPTCALLING.out.td_plot
}
