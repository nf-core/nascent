/*
 * Run parametertuning optionally, otherwise just run transcript calling
 */

include { GROHMM_PARAMETERTUNING } from '../../../modules/local/grohmm/parametertuning/main.nf'
include { GROHMM_TRANSCRIPTCALLING } from '../../../modules/local/grohmm/transcriptcalling/main.nf'

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

    // if(!tuning_file) {
        // Run transcriptcalling eval for each tuning param
        // Should avoid a tuning file with a row for everything
        // 5..45 by 5 for UTS is what we had currently
        ch_uts = channel.fromList((params.grohmm_min_uts..params.grohmm_max_uts).step(5))
        // -100..-400 by 50 for LtProbB
        ch_ltprobb = channel.fromList((params.grohmm_min_ltprobb..params.grohmm_max_ltprobb).step(50)).view()
        GROHMM_PARAMETERTUNING (
            bams_bais,
            gtf,
            ch_uts,
            ch_ltprobb,
        )
            .tuning
            .collectFile(
                name: "${params.outdir}/transcript_identification/grohmm/tuning.csv",
                keepHeader: true,
                skip: 1,
                newLine: true,
            )
            .set { tuning }

        ch_bams_bais_tuning = bams_bais.join(tuning, by: [0])

        ch_versions = ch_versions.mix(GROHMM_PARAMETERTUNING.out.versions.first())
    // } else {
        // If a tuning file is provided, run transcriptcalling once
        // NOTE This doesn't really handle multiple "groups well"
        // ch_bams_bais_tuning = bams_bais.join(tuning_file)
    // }

    GROHMM_TRANSCRIPTCALLING (
        ch_bams_bais_tuning,
        gtf,
    )
    ch_versions = ch_versions.mix(GROHMM_TRANSCRIPTCALLING.out.versions.first())

    emit:
    transcripts = GROHMM_TRANSCRIPTCALLING.out.transcripts
    bed = GROHMM_TRANSCRIPTCALLING.out.transcripts_bed
    td_plot = GROHMM_TRANSCRIPTCALLING.out.td_plot

    versions = ch_versions
}
