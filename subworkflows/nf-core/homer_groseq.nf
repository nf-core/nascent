/*
 * Identify transcripts with homer
 */

include { BEDTOOLS_BAMTOBED           } from '../../modules/nf-core/modules/bedtools/bamtobed/main'
include { HOMER_MAKETAGDIRECTORY      } from '../../modules/nf-core/modules/homer/maketagdirectory/main'
include { HOMER_MAKEUCSCFILE          } from '../../modules/nf-core/modules/homer/makeucscfile/main'
include { HOMER_FINDPEAKS             } from '../../modules/nf-core/modules/homer/findpeaks/main'

workflow HOMER_GROSEQ {
    take:
    bam // channel: [ val(meta), [ reads ] ]
    fasta //    file: /path/to/bwa/index/

    main:

    ch_versions = Channel.empty()

    // HACK Covert to bed files because homer doesn't have samtools included
    BEDTOOLS_BAMTOBED ( bam )
    BEDTOOLS_BAMTOBED.out.bed.map {
        meta, bed ->
        fmeta = meta.findAll { it.key != 'read_group' }
        fmeta.id = fmeta.id.split('_')[0..-2].join('_')
        [ fmeta, bed ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_grouped_bed }


    /*
    * Create a Tag Directory From The GRO-Seq experiment
    */
    HOMER_MAKETAGDIRECTORY ( ch_grouped_bed, fasta )
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
