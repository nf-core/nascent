include { PRESEQ_CCURVE } from '../../modules/nf-core/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP } from '../../modules/nf-core/preseq/lcextrap/main'
include { BBMAP_PILEUP } from '../../modules/nf-core/bbmap/pileup/main'

include { BAM_RSEQC  } from '../../subworkflows/nf-core/bam_rseqc'

workflow QUALITY_CONTROL {
    take:
    bam_bai
    bed

    main:

    bam = bam_bai.map{ [ it[0], it[1] ] }

    ch_versions = Channel.empty()

    PRESEQ_CCURVE ( bam )
    ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions.first())

    PRESEQ_LCEXTRAP ( bam )
    ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

    // TODO Set this in a param?
    rseqc_modules = ['read_duplication', 'read_distribution', 'infer_experiment']
    BAM_RSEQC ( bam_bai, bed, rseqc_modules )
    ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

    BBMAP_PILEUP ( bam )
    ch_versions = ch_versions.mix(BBMAP_PILEUP.out.versions.first())

    emit:
    preseq_ccurve = PRESEQ_CCURVE.out.c_curve
    preseq_lcextrap = PRESEQ_LCEXTRAP.out.lc_extrap

    inferexperiment_txt = BAM_RSEQC.out.inferexperiment_txt
    readdistribution_txt = BAM_RSEQC.out.readdistribution_txt
    readduplication_seq_xls = BAM_RSEQC.out.readduplication_seq_xls
    readduplication_pos_xls = BAM_RSEQC.out.readduplication_pos_xls
    readduplication_pdf = BAM_RSEQC.out.readduplication_pdf
    readduplication_rscript = BAM_RSEQC.out.readduplication_rscript

    pileup_stats = BBMAP_PILEUP.out.covstats
    pileup_hist = BBMAP_PILEUP.out.hist

    versions = ch_versions
}
