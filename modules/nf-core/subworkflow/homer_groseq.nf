/*
 * Identify transcripts with homer
 */

include { HOMER_MAKETAGDIRECTORY      } from '../software/homer/maketagdirectory/main' addParams( options: [:] )
include { HOMER_MAKEUCSCFILE          } from '../software/homer/makeucscfile/main'     addParams( options: [:] )
include { HOMER_FINDPEAKS             } from '../software/homer/findpeaks/main'        addParams( options: [:] )

workflow HOMER_GROSEQ {
    take:
    bam // channel: [ val(meta), [ reads ] ]
    fasta //    file: /path/to/bwa/index/

    main:
    /*
     * Create a Tag Directory From The GRO-Seq experiment
     */
    HOMER_MAKETAGDIRECTORY ( bam, fasta )

    emit:
    homer_version      = HOMER_MAKETAGDIRECTORY.out.version            //    path: *.version.txt
}
