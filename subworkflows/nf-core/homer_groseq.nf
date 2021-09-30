/*
 * Identify transcripts with homer
 */

params.maketagdirectory_options = [:]
params.makeucscfile_options     = [:]
params.findpeaks_options        = [args: '-style factor']

include { HOMER_MAKETAGDIRECTORY      } from '../../modules/nf-core/modules/homer/maketagdirectory/main' addParams( options: params.maketagdirectory_options )
include { HOMER_MAKEUCSCFILE          } from '../../modules/nf-core/modules/homer/makeucscfile/main'     addParams( options: params.makeucscfile_options )
include { HOMER_FINDPEAKS             } from '../../modules/nf-core/modules/homer/findpeaks/main'        addParams( options: params.findpeaks_options )

workflow HOMER_GROSEQ {
    take:
    bam // channel: [ val(meta), [ reads ] ]
    fasta //    file: /path/to/bwa/index/

    main:
    /*
     * Create a Tag Directory From The GRO-Seq experiment
     */
    HOMER_MAKETAGDIRECTORY ( bam, fasta )

    /*
     * Creating UCSC Visualization Files
     */
    HOMER_MAKEUCSCFILE ( HOMER_MAKETAGDIRECTORY.out.tagdir )

    /*
     * Find transcripts directly from GRO-Seq
     */
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )

    emit:
    homer_version      = HOMER_MAKETAGDIRECTORY.out.version //    path: *.version.txt
    tag_dir            = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), [ tag_dir ] ]

    bed_graph          = HOMER_MAKEUCSCFILE.out.bedGraph    // channel: [ val(meta), [ tag_dir/*ucsc.bedGraph.gz ] ]

    peaks              = HOMER_FINDPEAKS.out.txt            // channel: [ val(meta), [ *peaks.txt ] ]
}
