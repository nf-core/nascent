/*
 * Create bigWig and bedGraph files
 */

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PLUS
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MINUS } from '../../modules/nf-core/bedtools/genomecov/main'

include {
    DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_PLUS
    DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_MINUS } from '../../modules/nf-core/deeptools/bamcoverage/main'

workflow COVERAGE_GRAPHS {
    take:
    bam
    bai
    sizes
    fasta
    fai

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


    bam.join(bai, by: [0], remainder: true).set { ch_bam_bai }

    DEEPTOOLS_BAMCOVERAGE_PLUS (
        ch_bam_bai,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_PLUS.out.versions.first())

    DEEPTOOLS_BAMCOVERAGE_MINUS (
        ch_bam_bai,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_MINUS.out.versions.first())

    ch_plus_minus = DEEPTOOLS_BAMCOVERAGE_PLUS.out.bigwig.join(DEEPTOOLS_BAMCOVERAGE_MINUS.out.bigwig)

    emit:
    plus_bedGraph = BEDTOOLS_GENOMECOV_PLUS.out.genomecov
    minus_bedGraph = BEDTOOLS_GENOMECOV_MINUS.out.genomecov

    plus_minus = ch_plus_minus

    versions = ch_versions
}
