/*
 * Alignment with BWA
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { UNTAR             } from '../process/untar'                            addParams( options: params.index_options    )
include { BWA_INDEX         } from '../../nf-core/software/bwa/index/main'       addParams( options: params.index_options    )
include { BWA_MEM           } from '../../nf-core/software/bwa/mem/main'         addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS } from '../../nf-core/subworkflow/bam_sort_samtools' addParams( options: params.samtools_options )

workflow ALIGN_BWA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index //    file: /path/to/bwa/index/
    fasta //    file: /path/to/genome.fasta
    gtf   //    file: /path/to/genome.gtf

    main:
    /*
     * Uncompress BWA index or generate from scratch if required
    */
    if (index) {
        if (index.endsWith('.tar.gz')) {
            ch_index = UNTAR ( index ).untar
        } else {
            ch_index = file(index)
        }
    } else {
        ch_index = BWA_INDEX ( fasta, gtf ).index
    }

    /*
     * Map reads with BWA
     */
    BWA_MEM ( reads, ch_index, gtf )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( BWA_ALIGN.out.bam )

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
