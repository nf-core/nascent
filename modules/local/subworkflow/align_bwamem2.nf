/*
 * Alignment with BWA
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { UNTAR             } from '../process/untar'                            addParams( options: params.index_options    )
include { BWAMEM2_MEM       } from '../../nf-core/software/bwamem2/mem/main'         addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from '../../nf-core/subworkflow/bam_sort_samtools' addParams( options: params.samtools_options )

workflow ALIGN_BWAMEM2 {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/

    main:
    /*
     * Map reads with BWA
     */
    BWAMEM2_MEM ( reads, index )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( BWAMEM2_MEM.out.bam )

    emit:
    orig_bam         = BWAMEM2_MEM.out.bam                // channel: [ val(meta), bam ]
    bwa_version      = BWAMEM2_MEM.out.version            //    path: *.version.txt

    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
