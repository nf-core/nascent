//
// Alignment with BWAMEM2
//

include { BWAMEM2_MEM       } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { BAM_SORT_SAMTOOLS } from '../bam_sort_samtools/main'

workflow ALIGN_BWAMEM2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWAMEM2
    //
    BWAMEM2_MEM ( reads, index, true )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( BWAMEM2_MEM.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    orig_bam         = BWAMEM2_MEM.out.bam            // channel: [ val(meta), bam ]

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                      // channel: [ versions.yml ]
}
