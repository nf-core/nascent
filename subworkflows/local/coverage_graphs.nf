/*
 * TODO
 */

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PLUS
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MINUS } from '../../modules/nf-core/modules/bedtools/genomecov/main'

workflow COVERAGE_GRAPHS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    ch_genomecov_bam = bam.combine(Channel.from(1))

    BEDTOOLS_GENOMECOV_PLUS (
        ch_genomecov_bam,
        [],
        'bedGraph'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_PLUS.out.versions.first())

    BEDTOOLS_GENOMECOV_MINUS (
        ch_genomecov_bam,
        [],
        'bedGraph'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_MINUS.out.versions.first())

    ch_plus_minus = BEDTOOLS_GENOMECOV_PLUS.out.genomecov.join(BEDTOOLS_GENOMECOV_MINUS.out.genomecov)

    emit:
    plus_bedGraph = BEDTOOLS_GENOMECOV_PLUS.out.genomecov
    minus_bedGraph = BEDTOOLS_GENOMECOV_MINUS.out.genomecov

    plus_minus = ch_plus_minus

    versions = ch_versions                      // channel: [ versions.yml ]
}
