/*
 * Alignment with BWA
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { BWA_MEM           } from '../../modules/nf-core/modules/bwa/mem/main'   addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from './bam_sort_samtools'  addParams( options: params.samtools_options )

workflow ALIGN_BWA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/

    main:
    /*
     * Map reads with BWA
     */
    BWA_MEM ( reads, index )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( BWA_MEM.out.bam )

    emit:
    orig_bam         = BWA_MEM.out.bam                // channel: [ val(meta), bam ]
    bwa_version      = BWA_MEM.out.version            //    path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
