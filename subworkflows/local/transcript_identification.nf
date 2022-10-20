include { GROHMM } from './grohmm'

include { PINTS_CALLER } from '../../modules/nf-core/pints/caller/main'
include { CAT_CAT } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_MERGE } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_FILTER } from '../../modules/nf-core/bedtools/intersect/main'

include { HOMER_GROSEQ } from '../nf-core/homer/groseq/main'

workflow TRANSCRIPT_INDENTIFICATION {
    take:
    group_bams
    gtf
    fasta

    main:

    ch_versions = Channel.empty()
    ch_identification_bed = Channel.empty()

    ch_tuning_file = params.tuning_file ? file(params.tuning_file, checkIfExists: true) : file("${projectDir}/assets/tuningparamstotest.csv")
    GROHMM ( group_bams, gtf, ch_tuning_file )
    ch_identification_bed = ch_identification_bed.mix(GROHMM.out.bed)
    ch_versions = ch_versions.mix(GROHMM.out.versions.first())


    HOMER_GROSEQ ( group_bams, fasta )
    ch_identification_bed = ch_identification_bed.mix(HOMER_GROSEQ.out.bed)
    ch_versions = ch_versions.mix(HOMER_GROSEQ.out.versions.first())


    // TODO https://github.com/hyulab/PINTS/issues/15
    PINTS_CALLER ( group_bams )
    ch_versions = ch_versions.mix(PINTS_CALLER.out.versions.first())
    // HACK Not sure if this is as good as reporting all of them, but it should
    // reduce the overall noise.
    CAT_CAT ( PINTS_CALLER.out.bidirectional_TREs )
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())
    BEDTOOLS_SORT ( CAT_CAT.out.file_out, "bed" )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions.first())
    BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )
    ch_identification_bed = ch_identification_bed.mix(BEDTOOLS_MERGE.out.bed)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    ch_filter_bed = Channel.from(params.filter_bed)
    BEDTOOLS_INTERSECT_FILTER ( ch_identification_bed.combine(ch_filter_bed), "bed" )
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_FILTER.out.versions.first())

    // Use non-filtered bed files if we skip filtering
    if(params.filter_bed) {
        ch_identification_bed = BEDTOOLS_INTERSECT_FILTER.out.intersect
    }

    ch_identification_bed
        // Drop any empty bed files
        .filter { meta, bed -> bed.size() > 0 }
        .set { ch_identification_bed_clean }

    emit:
    grohmm_td_plot = GROHMM.out.td_plot
    homer_peaks = HOMER_GROSEQ.out.peaks
    homer_tagdir = HOMER_GROSEQ.out.tagdir

    transcript_beds = ch_identification_bed_clean

    versions = ch_versions
}
