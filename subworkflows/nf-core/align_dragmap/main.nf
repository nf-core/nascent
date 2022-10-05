//
// Alignment with dragmap
//

include { DRAGMAP_ALIGN } from '../../../modules/nf-core/dragmap/align/main'
include { BAM_SORT_SAMTOOLS } from '../bam_sort_samtools/main'

workflow ALIGN_DRAGMAP {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWA MEM
    //
    DRAGMAP_ALIGN ( reads, index, true )
    ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( DRAGMAP_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    orig_bam         = DRAGMAP_ALIGN.out.bam                // channel: [ val(meta), bam ]

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                      // channel: [ versions.yml ]
}
