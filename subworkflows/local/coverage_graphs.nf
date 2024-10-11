/*
 * Create bigWig and bedGraph files
 */

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_PLUS
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MINUS } from '../../modules/nf-core/bedtools/genomecov/main'

include {
    DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_PLUS
    DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_MINUS } from '../../modules/nf-core/deeptools/bamcoverage/main'

include { DREG_PREP } from '../../modules/local/dreg_prep/main'

workflow COVERAGE_GRAPHS {
    take:
    bam_bai
    sizes
    fasta
    fai

    main:

    bam = bam_bai.map{ [ it[0], it[1] ] }

    ch_genomecov_bam = bam.combine(Channel.from(1))

    BEDTOOLS_GENOMECOV_PLUS (
        ch_genomecov_bam,
        sizes,
        'bedGraph',
        true
    )

    BEDTOOLS_GENOMECOV_MINUS (
        ch_genomecov_bam,
        sizes,
        'bedGraph',
        true
    )

    DEEPTOOLS_BAMCOVERAGE_PLUS (
        bam_bai,
        fasta,
        fai
    )

    DEEPTOOLS_BAMCOVERAGE_MINUS (
        bam_bai,
        fasta,
        fai
    )

    ch_plus_minus = DEEPTOOLS_BAMCOVERAGE_PLUS.out.bigwig.join(DEEPTOOLS_BAMCOVERAGE_MINUS.out.bigwig)

    DREG_PREP (
        bam_bai,
        sizes
    )

    emit:
    plus_bedGraph = BEDTOOLS_GENOMECOV_PLUS.out.genomecov
    minus_bedGraph = BEDTOOLS_GENOMECOV_MINUS.out.genomecov

    plus_minus = ch_plus_minus
}
