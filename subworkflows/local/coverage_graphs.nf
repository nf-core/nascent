/*
 * TODO
 */

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PLUS
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MINUS } from '../../modules/nf-core/modules/bedtools/genomecov/main'

include {
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_PLUS
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_MINUS } from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'

workflow COVERAGE_GRAPHS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    sizes

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

    UCSC_BEDGRAPHTOBIGWIG_PLUS (
        BEDTOOLS_GENOMECOV_PLUS.out.genomecov,
        sizes
    )

    UCSC_BEDGRAPHTOBIGWIG_MINUS (
        BEDTOOLS_GENOMECOV_MINUS.out.genomecov,
        sizes
    )

    ch_plus_minus = UCSC_BEDGRAPHTOBIGWIG_PLUS.out.bigwig.join(UCSC_BEDGRAPHTOBIGWIG_MINUS.out.bigwig)

    emit:
    plus_bedGraph = BEDTOOLS_GENOMECOV_PLUS.out.genomecov
    minus_bedGraph = BEDTOOLS_GENOMECOV_MINUS.out.genomecov

    plus_minus = ch_plus_minus

    versions = ch_versions                      // channel: [ versions.yml ]
}
