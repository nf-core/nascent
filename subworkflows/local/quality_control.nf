include { PRESEQ_CCURVE } from '../../modules/nf-core/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP } from '../../modules/nf-core/preseq/lcextrap/main'
include { RSEQC_READDISTRIBUTION } from '../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_READDUPLICATION } from '../../modules/nf-core/rseqc/readduplication/main'
include { RSEQC_INFEREXPERIMENT } from '../../modules/nf-core/rseqc/inferexperiment/main'
include { BBMAP_PILEUP } from '../../modules/nf-core/bbmap/pileup/main'

workflow QUALITY_CONTROL {
    take:
    bam
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

    RSEQC_READDUPLICATION ( bam )
    ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())

    RSEQC_INFEREXPERIMENT (
        bam,
        bed
    )
    ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())

    BBMAP_PILEUP ( bam )
    ch_versions = ch_versions.mix(BBMAP_PILEUP.out.versions.first())

    emit:
    preseq_ccurve = PRESEQ_CCURVE.out.c_curve
    preseq_lcextrap = PRESEQ_LCEXTRAP.out.lc_extrap

    readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
    readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
    readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
    readduplication_pdf = RSEQC_READDUPLICATION.out.pdf
    readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
    inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.txt

    pileup_stats = BBMAP_PILEUP.out.covstats
    pileup_hist = BBMAP_PILEUP.out.hist

    versions = ch_versions
}
