include { GROHMM                } from './grohmm'

include { PINTS_CALLER                                            } from '../../modules/nf-core/pints/caller/main'
include { BEDTOOLS_MERGE                                          } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_FILTER         } from '../../modules/nf-core/bedtools/intersect/main'

include { HOMER_GROSEQ          } from '../nf-core/homer/groseq/main'

workflow TRANSCRIPT_INDENTIFICATION {
    take:
    group_bams // channel: [ val(meta), [ bams ] ]
    gtf
    fasta

    main:

    ch_identification_bed = Channel.empty()

    ch_grohmm_multiqc = Channel.empty()
    ch_tuning_file = params.tuning_file ? file(params.tuning_file, checkIfExists: true) : file("${projectDir}/assets/tuningparamstotest.csv")
    GROHMM (
        group_bams,
        gtf,
        ch_tuning_file
    )
    ch_grohmm_multiqc = GROHMM.out.td_plot.collect()
    ch_identification_bed = ch_identification_bed.mix(GROHMM.out.bed)


    ch_homer_multiqc = Channel.empty()
    HOMER_GROSEQ (
        group_bams,
        fasta
    )
    ch_homer_multiqc = HOMER_GROSEQ.out.peaks
    ch_homer_multiqc = ch_homer_multiqc.mix(HOMER_GROSEQ.out.tagdir)
    ch_identification_bed = ch_identification_bed.mix(HOMER_GROSEQ.out.bed)


    // TODO Merge technical replicates
    PINTS_CALLER (
        group_bams
    )
    // HACK Not sure if this is as good as reporting all of them, but it should
    // reduce the overall noise.
    BEDTOOLS_MERGE (
        PINTS_CALLER.out.bidirectional_TREs
    )
    ch_identification_bed = ch_identification_bed.mix(BEDTOOLS_MERGE.out.bed)

    // TODO Support gzipped bed files
    ch_filter_bed = Channel.from(params.filter_bed)
    BEDTOOLS_INTERSECT_FILTER (
        ch_identification_bed.combine(ch_filter_bed),
        "bed"
    )

    // Use non-filtered bed files if we skip filtering
    if(params.filter_bed) {
        ch_identification_bed = BEDTOOLS_INTERSECT_FILTER.out.intersect
    }

    ch_identification_bed
        // Drop any empty bed files
        .filter { meta, bed -> bed.size() > 0 }
        .set { ch_identification_bed_clean }

    // Gather versions of all tools used
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(GROHMM.out.versions.first())
    ch_versions = ch_versions.mix(HOMER_GROSEQ.out.versions.first())
    ch_versions = ch_versions.mix(PINTS_CALLER.out.versions.first())
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_FILTER.out.versions.first())

    emit:
    grohmm_td_plot = GROHMM.out.td_plot
    homer_peaks = HOMER_GROSEQ.out.peaks
    homer_tagdir = HOMER_GROSEQ.out.tagdir

    transcript_beds = ch_identification_bed_clean

    versions = ch_versions
}
