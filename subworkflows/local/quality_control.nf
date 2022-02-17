/*
 * TODO
 */

include { PRESEQ_LCEXTRAP } from '../../modules/nf-core/modules/preseq/lcextrap/main'

workflow QUALITY_CONTROL {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    PRESEQ_LCEXTRAP ( bam )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    emit:
    preseq_lcextrap = PRESEQ_LCEXTRAP.out.ccurve

    versions        = ch_versions                      // channel: [ versions.yml ]
}
