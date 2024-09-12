/*
 * Run parametertuning optionally, otherwise just run transcript calling
 */

include { GROHMM_TRANSCRIPTCALLING } from '../../../modules/local/grohmm/transcriptcalling/main.nf'
include { GROHMM_PARAMETERTUNING } from '../../../modules/local/grohmm/parametertuning/main.nf'

/*
 * Note meta refers to all merged files
 */
workflow GROHMM {
    take:
    bams_bais
    gtf
    tuning_file

    main:

    ch_versions = Channel.empty()

    ch_tuning = Channel.empty()

    if(!params.skip_tuning) {
        GROHMM_PARAMETERTUNING (
            bams_bais,
            gtf,
            tuning_file
        )
        ch_tuning = GROHMM_PARAMETERTUNING.out.tuning
        ch_bams_bais_tuning = bams_bais.join(ch_tuning)
        ch_versions = ch_versions.mix(GROHMM_PARAMETERTUNING.out.versions.first())
    } else {
        ch_bams_bais_tuning = bams_bais.join(ch_tuning)
    }

    GROHMM_TRANSCRIPTCALLING (
        [ch_bams_bais_tuning, []],
        gtf,
    )
    ch_versions = ch_versions.mix(GROHMM_TRANSCRIPTCALLING.out.versions.first())

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    td_plot = GROHMM_TRANSCRIPTCALLING.out.td_plot

    versions = ch_versions
}
