/*
 * TODO
 */

include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/modules/bedtools/genomecov/main'

workflow COVERAGE_GRAPHS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    ch_genomecov_bam = bam.combine(Channel.from(1))

    BEDTOOLS_GENOMECOV (
        ch_genomecov_bam,
        [],
        'bedGraph'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    emit:
    bedGraph = BEDTOOLS_GENOMECOV.out.genomecov

    versions = ch_versions                      // channel: [ versions.yml ]
}
