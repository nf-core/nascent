/*
 * Calls Transcripts and Transcript Start Sites and various cleaning steps
 */

include { GROHMM } from './grohmm/main'
include { HOMER_GROSEQ } from '../nf-core/homer/groseq/main'
include { PINTS_CALLER } from '../../modules/nf-core/pints/caller/main'
include { CAT_CAT } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_MERGE } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_FILTER } from '../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT } from '../../modules/nf-core/bedtools/intersect/main'

workflow TRANSCRIPT_INDENTIFICATION {
    take:
    group_bam_bai
    gtf
    fasta
    chrom_sizes

    main:

    ch_identification_bed = Channel.empty()

    grohmm_td_plot = Channel.empty()
    if(!params.skip_grohmm && params.assay_type == "GROseq") {
        GROHMM ( group_bam_bai, gtf )
        ch_identification_bed = ch_identification_bed.mix(GROHMM.out.bed)
        grohmm_td_plot = GROHMM.out.td_plot
    }


    homer_peaks = Channel.empty()
    homer_tagdir = Channel.empty()
    if(params.assay_type == "GROseq") {
        group_bam = group_bam_bai.map { meta, bam, bai -> [meta, bam] }
        HOMER_GROSEQ ( group_bam, fasta )
        ch_identification_bed = ch_identification_bed.mix(HOMER_GROSEQ.out.bed)
        homer_peaks = HOMER_GROSEQ.out.peaks
        homer_tagdir = HOMER_GROSEQ.out.tagdir
    }


    // TODO https://github.com/hyulab/PINTS/issues/15
    PINTS_CALLER ( group_bam_bai, params.assay_type )
    // HACK Not sure if this is as good as reporting all of them, but it should
    // reduce the overall noise.
    CAT_CAT ( PINTS_CALLER.out.bidirectional_TREs )
    BEDTOOLS_SORT ( CAT_CAT.out.file_out, [] )
    BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )
    ch_identification_bed = ch_identification_bed.mix(BEDTOOLS_MERGE.out.bed)

    if(params.filter_bed) {
        ch_filter_bed = Channel.from(params.filter_bed)
        BEDTOOLS_INTERSECT_FILTER ( ch_identification_bed.combine(ch_filter_bed), chrom_sizes.map { [ [:], it ] } )
        ch_identification_bed = BEDTOOLS_INTERSECT_FILTER.out.intersect
    }
    if(params.intersect_bed) {
        ch_intersect_bed = Channel.from(params.intersect_bed)
        BEDTOOLS_INTERSECT ( ch_identification_bed.combine(ch_intersect_bed), chrom_sizes.map { [ [:], it ] } )
        ch_identification_bed = BEDTOOLS_INTERSECT.out.intersect
    }

    ch_identification_bed
        // Drop any empty bed files
        .filter { meta, bed -> bed.size() > 0 }
        .set { ch_identification_bed_clean }

    emit:
    grohmm_td_plot
    homer_peaks
    homer_tagdir

    transcript_beds = ch_identification_bed_clean
}
