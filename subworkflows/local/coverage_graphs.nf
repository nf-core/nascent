/*
 * Create bigWig and bedGraph files
 */

include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_PLUS  } from '../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_MINUS } from '../../modules/nf-core/deeptools/bamcoverage/main'

include { DREG_PREP                                            } from '../../modules/local/dreg_prep/main'

workflow COVERAGE_GRAPHS {
    take:
    bam_bai
    sizes
    fasta
    fai

    main:

    bam = bam_bai.map { [it[0], it[1]] }

    ch_versions = Channel.empty()

    ch_genomecov_bam = bam.combine(Channel.from(1))

    DEEPTOOLS_BAMCOVERAGE_PLUS(
        bam_bai,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_PLUS.out.versions.first())

    DEEPTOOLS_BAMCOVERAGE_MINUS(
        bam_bai,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_MINUS.out.versions.first())

    ch_plus_minus = DEEPTOOLS_BAMCOVERAGE_PLUS.out.bigwig.join(DEEPTOOLS_BAMCOVERAGE_MINUS.out.bigwig)

    DREG_PREP(
        bam_bai,
        sizes,
        params.assay_type
    )

    emit:
    plus_minus     = ch_plus_minus
    versions       = ch_versions
}
