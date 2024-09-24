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
    // TODO Support rerunning with a tuning file
    // tuning_file

    main:

    ch_versions = Channel.empty()

    // Run transcriptcalling eval for each tuning param
    // Should avoid a tuning file with a row for everything
    // 5..45 by 5 for UTS is what we had currently
    ch_uts = channel.fromList((params.grohmm_min_uts..params.grohmm_max_uts).step(5))
    // -100..-400 by 50 for LtProbB
    ch_ltprobb = channel.fromList((params.grohmm_min_ltprobb..params.grohmm_max_ltprobb).step(50))
    ch_bams_bais_uts_ltprobb =
        bams_bais
        .combine(ch_uts)
        .combine(ch_ltprobb)

    GROHMM_PARAMETERTUNING (
        ch_bams_bais_uts_ltprobb,
        gtf,
    )
        .tuning
        .collectFile(
            name: "${params.outdir}/transcript_identification/grohmm/tuning.csv",
            keepHeader: true,
            skip: 1,
            newLine: true,
        )
        .set { ch_tuning }
    ch_versions = ch_versions.mix(GROHMM_PARAMETERTUNING.out.versions.first())

    ch_bams_bais_tuning = bams_bais.join(ch_tuning, by: [0])


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
