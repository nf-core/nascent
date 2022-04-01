//
// Alignment with BWA MEM
//

include { BWA_MEM           } from '../../../modules/nf-core/modules/bwa/mem/main'
include { BAM_SORT_SAMTOOLS } from '../bam_sort_samtools/main'

workflow ALIGN_BWA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with BWA MEM
    //
    BWA_MEM ( reads, index, true )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_SAMTOOLS ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    orig_bam         = BWA_MEM.out.bam                // channel: [ val(meta), bam ]

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                      // channel: [ versions.yml ]
}
