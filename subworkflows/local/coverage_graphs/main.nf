/*
 * Create bigWig and bedGraph files
 */

include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_PLUS  } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_MINUS } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { FLIP_STRAND                                          } from '../../../modules/local/flip_strand'
include { PINTS_VISUALIZER                                     } from '../../../modules/local/pints/visualizer'

workflow COVERAGE_GRAPHS {
    take:
    bam_bai
    sizes
    fasta
    fai

    main:

    bam = bam_bai.map { [it[0], it[1]] }

    ch_versions = Channel.empty()

    DEEPTOOLS_BAMCOVERAGE_PLUS(
        bam_bai,
        fasta,
        fai,
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_PLUS.out.versions.first())

    DEEPTOOLS_BAMCOVERAGE_MINUS(
        bam_bai,
        fasta,
        fai,
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_MINUS.out.versions.first())


    ch_plus_minus = DEEPTOOLS_BAMCOVERAGE_PLUS.out.bigwig.join(DEEPTOOLS_BAMCOVERAGE_MINUS.out.bigwig)

    if (params.assay_type in ["PROseq", "PROcap", "CoPRO"]) {
        FLIP_STRAND(
            ch_plus_minus,
        )
        ch_versions = ch_versions.mix(FLIP_STRAND.out.versions.first())
        ch_plus_minus = FLIP_STRAND.out.flipped_bigwig
    }

    emit:
    plus_minus = ch_plus_minus
    versions   = ch_versions
}
