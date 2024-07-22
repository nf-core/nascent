//
// Alignment with dragmap
//

include { DRAGMAP_ALIGN } from '../../../modules/nf-core/dragmap/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_DRAGMAP {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam    // boolean (mandatory): true or false
    ch_fasta        // channel (optional) : [ val(meta3), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with dragmap
    //

    DRAGMAP_ALIGN ( ch_reads, ch_index, ch_fasta, val_sort_bam )
    ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //

    BAM_SORT_STATS_SAMTOOLS ( DRAGMAP_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = DRAGMAP_ALIGN.out.bam                // channel: [ val(meta), path(bam) ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
