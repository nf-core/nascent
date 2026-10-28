//
// Alignment with MiniBWA
//

include { MINIBWA_MAP             } from '../../../modules/local/minibwa/map/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_MINIBWA {
    take:
    ch_reads // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam // boolean (mandatory): true or false
    ch_fasta // channel (optional) : [ val(meta3), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    MINIBWA_MAP(ch_reads, ch_index, ch_fasta, val_sort_bam)
    ch_versions = ch_versions.mix(MINIBWA_MAP.out.versions.first())

    BAM_SORT_STATS_SAMTOOLS(MINIBWA_MAP.out.bam, ch_fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = MINIBWA_MAP.out.bam
    log_out  = MINIBWA_MAP.out.log
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    versions = ch_versions
}
