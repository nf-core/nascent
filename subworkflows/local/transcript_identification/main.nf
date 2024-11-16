/*
 * Calls Transcripts and Transcript Start Sites and various cleaning steps
 */

include { GROHMM                                          } from '../grohmm'
include { HOMER_GROSEQ                                    } from '../../nf-core/homer/groseq/main'
include { PINTS_CALLER                                    } from '../../../modules/nf-core/pints/caller/main'
include { SAMTOOLS_MERGE                                  } from '../../../modules/nf-core/samtools/merge/main'
include { CAT_CAT                                         } from '../../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_MERGE                                  } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                                   } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_FILTER } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT                              } from '../../../modules/nf-core/bedtools/intersect/main'

workflow TRANSCRIPT_INDENTIFICATION {
    take:
    group_bam_bai
    gxf
    fasta
    chrom_sizes
    uniqmap

    main:

    ch_versions = Channel.empty()
    ch_identification_bed = Channel.empty()

    grohmm_td_plot = Channel.empty()
    if (!params.skip_grohmm && params.assay_type == "GROseq") {
        GROHMM(group_bam_bai, gxf)
        ch_identification_bed = ch_identification_bed.mix(GROHMM.out.bed)
        grohmm_td_plot = GROHMM.out.td_plot
        ch_versions = ch_versions.mix(GROHMM.out.versions.first())
    }


    homer_peaks = Channel.empty()
    homer_tagdir = Channel.empty()
    if (params.assay_type == "GROseq") {
        group_bam = group_bam_bai.map { meta, bam, bai -> [meta, bam] }
        HOMER_GROSEQ(group_bam, fasta, uniqmap)
        ch_identification_bed = ch_identification_bed.mix(HOMER_GROSEQ.out.bed)
        homer_peaks = HOMER_GROSEQ.out.peaks
        homer_tagdir = HOMER_GROSEQ.out.tagdir
        ch_versions = ch_versions.mix(HOMER_GROSEQ.out.versions.first())
    }

    // Scatter the chromosomes
    // TODO Might turn this into a param
    skip_chr = ["chrY", "_random", "chrUn_", "chrEBV", "chrM"]
    ch_chr = fasta
        .splitFasta(record: [id: true])
        .map { record -> record.id }
        .filter { !(it in skip_chr) }

    // NOTE https://github.com/hyulab/PINTS/issues/15
    // We want to run PINTS once per biological/technical sample.
    // Let Nextflow handle the parallelization Per the docs
    // Per the docs:
    // {SID} will be replaced with the number of samples that peaks are called from,
    // if you only provide PINTS with one sample, then {SID} will be replaced with 1,
    // if you try to use PINTS with three replicates (--bam-file A.bam B.bam C.bam),
    // then {SID} for peaks identified from A.bam will be replaced with 1.
    SAMTOOLS_MERGE(
        group_bam_bai.map { meta, bams, _bais -> [meta, bams] },
        [[], []],
        [[], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    PINTS_CALLER(
        SAMTOOLS_MERGE.out.bam.combine(ch_chr),
        params.assay_type
    )
    ch_versions = ch_versions.mix(PINTS_CALLER.out.versions.first())

    // Gather the chromosomes
    // TODO Tests don't seem to hit this because there's no bidirectional_TREs
    // Need to collect all of the beds for each chromosome/sample and concatenate them
    // Nextflow makes this super easy
    def ch_bidirectional_tres = PINTS_CALLER.out.unidirectional_TREs.groupTuple(by: [0]).map { meta, beds ->
        [meta, beds.flatten()]
    }

    // HACK Not sure if this is as good as reporting all of them, but it should
    // reduce the overall noise.
    CAT_CAT(ch_bidirectional_tres)
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())
    BEDTOOLS_SORT(CAT_CAT.out.file_out, [])
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions.first())
    BEDTOOLS_MERGE(BEDTOOLS_SORT.out.sorted)
    ch_identification_bed = ch_identification_bed.mix(BEDTOOLS_MERGE.out.bed)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    if (params.filter_bed) {
        ch_filter_bed = Channel.from(params.filter_bed)
        BEDTOOLS_INTERSECT_FILTER(ch_identification_bed.combine(ch_filter_bed.first()), chrom_sizes.map { [[:], it] })
        ch_identification_bed = BEDTOOLS_INTERSECT_FILTER.out.intersect
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_FILTER.out.versions.first())
    }
    if (params.intersect_bed) {
        ch_intersect_bed = Channel.from(params.intersect_bed)
        BEDTOOLS_INTERSECT(ch_identification_bed.combine(ch_intersect_bed.first()), chrom_sizes.map { [[:], it] })
        ch_identification_bed = BEDTOOLS_INTERSECT.out.intersect
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())
    }

    ch_identification_bed
        .filter { meta, bed -> bed.size() > 0 }
        .set { ch_identification_bed_clean }

    emit:
    grohmm_td_plot
    homer_peaks
    homer_tagdir
    transcript_beds = ch_identification_bed_clean
    versions        = ch_versions
}
