/*
 * TODO
 */

include { PRESEQ_CCURVE } from '../../modules/local/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP } from '../../modules/nf-core/modules/preseq/lcextrap/main'
include { RSEQC_READDISTRIBUTION } from '../../modules/nf-core/modules/rseqc/readdistribution/main'

workflow QUALITY_CONTROL {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    bed

    main:

    ch_versions = Channel.empty()

    PRESEQ_CCURVE ( bam )
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions.first())

    PRESEQ_LCEXTRAP ( bam )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    RSEQC_READDISTRIBUTION (
        bam,
        bed
    )
    ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())

    emit:
    preseq_ccurve = PRESEQ_CCURVE.out.ccurve
    preseq_lcextrap = PRESEQ_LCEXTRAP.out.ccurve

    readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
    versions        = ch_versions                      // channel: [ versions.yml ]
}
