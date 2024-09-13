/*
 * Run parametertuning optionally, otherwise just run transcript calling
 */

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

    // If a tuning file is provided, run transcriptcalling once
    if(tuning_file) {
        // TODO Find minimum
        // uts <- tune[which.min(tune$errorRate), "UTS"]
        // lt_probb <- tune[which.min(tune$errorRate), "LtProbB"]
        GROHMM_TRANSCRIPTCALLING(
            bam_bais,
            gtf,
            minimum_uts,
            minimum_ltprobb,
        )
    } else {
        // Run transcriptcalling eval for each tuning param
        // Should avoid a tuning file with a row for everything
        // 5..45 by 5 for UTS is what we had currently
        // -100..-400 by 50 for LtProbB
        GROHMM_PARAMETERTUNING (
            bam_bais,
            gtf,
            ch_uts,
            ch_ltprobb,
        )
        // TODO CollectFile the tuning
        // TODO Find the minimum values
        // TODO Need to decide if windowAnalysis is important
        // https://github.com/Functional-Genomics-Lab/groseq-analysis/blob/9b69519c41232fd653a2b2726e32d91b49abeb7e/research/groHMM2.R#L62C7-L62C21
        // If it is need to rerun transcriptcalling without it
    }



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
