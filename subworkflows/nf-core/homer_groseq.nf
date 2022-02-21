/*
 * Identify transcripts with homer
 */

include { HOMER_MAKETAGDIRECTORY      } from '../../modules/nf-core/modules/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE          } from '../../modules/nf-core/modules/homer/makeucscfile/main'
include { HOMER_FINDPEAKS             } from '../../modules/nf-core/modules/homer/findpeaks/main'

workflow HOMER_GROSEQ {
    take:
    bam // channel: [ val(meta), [ reads ] ]
    fasta //    file: /path/to/bwa/index/

    main:

    ch_versions = Channel.empty()

    /*
    * Create a Tag Directory From The GRO-Seq experiment
    */
    HOMER_MAKETAGDIRECTORY ( bam, fasta )
    ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions.first())

    /*
    * Creating UCSC Visualization Files
    */
    HOMER_MAKEUCSCFILE ( HOMER_MAKETAGDIRECTORY.out.tagdir )
    ch_versions = ch_versions.mix(HOMER_MAKEUCSCFILE.out.versions.first())

    /*
    * Find transcripts directly from GRO-Seq
    */
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )
    ch_versions = ch_versions.mix(HOMER_FINDPEAKS.out.versions.first())

    emit:
    tag_dir            = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), [ tag_dir ] ]
    bed_graph          = HOMER_MAKEUCSCFILE.out.bedGraph    // channel: [ val(meta), [ tag_dir/*ucsc.bedGraph.gz ] ]
    peaks              = HOMER_FINDPEAKS.out.txt            // channel: [ val(meta), [ *peaks.txt ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}
