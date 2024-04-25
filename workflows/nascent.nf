/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BED2SAF } from '../modules/local/bed2saf'

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_BWAMEM2 } from '../subworkflows/local/align_bwamem2/main'
include { ALIGN_DRAGMAP } from '../subworkflows/local/align_dragmap/main'
include { QUALITY_CONTROL } from '../subworkflows/local/quality_control.nf'
include { COVERAGE_GRAPHS } from '../subworkflows/local/coverage_graphs.nf'
include { TRANSCRIPT_INDENTIFICATION } from '../subworkflows/local/transcript_identification.nf'

include { FASTP } from '../modules/nf-core/fastp/main'
include {
    UNTAR as UNTAR_HISAT2_INDEX
    UNTAR as UNTAR_STAR_INDEX
} from '../modules/nf-core/untar/main'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'
include {
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_PREDICTED } from '../modules/nf-core/subread/featurecounts/main'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nascent_pipeline'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_ALIGN_BWA } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { FASTQ_ALIGN_BOWTIE2 } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { FASTQ_ALIGN_HISAT2 } from '../subworkflows/nf-core/fastq_align_hisat2/main'
include { FASTQ_ALIGN_STAR } from '../subworkflows/nf-core/fastq_align_star/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NASCENT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_gtf
    ch_gff
    ch_gene_bed
    ch_bwa_index
    ch_bwamem2_index
    ch_dragmap
    ch_bowtie2_index
    ch_hisat2_index
    ch_star_index

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //

    // HACK Rework this because of nf-validation
    def prepareToolIndices = []
    if (!params.skip_alignment) { prepareToolIndices << params.aligner }
    PREPARE_GENOME (
        prepareToolIndices,
        ch_fasta,
        ch_gtf,
        ch_gff,
        ch_gene_bed,
        ch_bwa_index,
        ch_bwamem2_index,
        ch_dragmap,
        ch_bowtie2_index,
        ch_hisat2_index,
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.first())
    ch_fasta = PREPARE_GENOME.out.fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ch_reads = Channel.empty()
    if(!params.skip_trimming) {
        FASTP ( ch_samplesheet, [], false, false )
        ch_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    } else {
        ch_reads = ch_samplesheet
    }

    //
    // SUBWORKFLOW: Alignment with BWA
    //
    ch_genome_bam = Channel.empty()
    ch_genome_bai = Channel.empty()
    ch_samtools_stats = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_star_multiqc = Channel.empty()
    ch_aligner_pca_multiqc = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    ch_bowtie2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bwa') {
        FASTQ_ALIGN_BWA (
            ch_reads,
            PREPARE_GENOME.out.bwa_index,
            false,
            ch_fasta,
        )
        ch_genome_bam = FASTQ_ALIGN_BWA.out.bam
        ch_genome_bai = FASTQ_ALIGN_BWA.out.bai
        ch_samtools_stats = FASTQ_ALIGN_BWA.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BWA.out.idxstats

        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions.first())
    } else if (!params.skip_alignment && params.aligner == 'bwamem2') {
        ALIGN_BWAMEM2 (
            ch_reads,
            PREPARE_GENOME.out.bwa_index,
            false,
            ch_fasta,
        )
        ch_genome_bam = ALIGN_BWAMEM2.out.bam
        ch_genome_bai = ALIGN_BWAMEM2.out.bai
        ch_samtools_stats = ALIGN_BWAMEM2.out.stats
        ch_samtools_flagstat = ALIGN_BWAMEM2.out.flagstat
        ch_samtools_idxstats = ALIGN_BWAMEM2.out.idxstats

        ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions)
    } else if (!params.skip_alignment && params.aligner == 'dragmap') {
        ALIGN_DRAGMAP (
            ch_reads,
            PREPARE_GENOME.out.dragmap,
            false,
            ch_fasta,
        )
        ch_genome_bam = ALIGN_DRAGMAP.out.bam
        ch_genome_bai = ALIGN_DRAGMAP.out.bai
        ch_samtools_stats = ALIGN_DRAGMAP.out.stats
        ch_samtools_flagstat = ALIGN_DRAGMAP.out.flagstat
        ch_samtools_idxstats = ALIGN_DRAGMAP.out.idxstats

        ch_versions = ch_versions.mix(ALIGN_DRAGMAP.out.versions)
    } else if (!params.skip_alignment && params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            ch_reads,
            PREPARE_GENOME.out.bowtie2_index,
            false,
            false,
            ch_fasta,
        )
        ch_genome_bam = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_genome_bai = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_samtools_stats = FASTQ_ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats

        ch_bowtie2_multiqc = FASTQ_ALIGN_BOWTIE2.out.log_out
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)
    } else if (!params.skip_alignment && params.aligner == 'hisat2') {
        if (ch_hisat2_index.endsWith('.tar.gz')) {
            ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], ch_hisat2_index ] ).untar
            ch_versions = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
        } else {
            // TODO Give the meta from basename or genome?
            ch_hisat2_index = [ [meta: "Genome"], file(ch_hisat2_index) ]
        }

        FASTQ_ALIGN_HISAT2 (
            ch_reads,
            ch_hisat2_index,
            [[:],[]],
            ch_fasta,
        )
        ch_genome_bam = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bai = FASTQ_ALIGN_HISAT2.out.bai
        ch_samtools_stats = FASTQ_ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats

        ch_HISAT2_multiqc = FASTQ_ALIGN_HISAT2.out.summary
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)
    } else if (!params.skip_alignment && params.aligner == 'star') {
        if(!ch_star_index) {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, [ [:], ch_gtf ] ).index
        } else if (ch_star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( [ [:], ch_star_index ] ).untar
            ch_versions = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            // TODO Give the meta from basename or genome?
            ch_star_index = [ [meta: "Genome"], file(ch_star_index) ]
        }

        FASTQ_ALIGN_STAR (
            ch_reads,
            ch_star_index,
            PREPARE_GENOME.out.gtf,
            false,
            '',
            '', // TODO params.seq_center ?:
            ch_fasta,
            []
        )
        ch_genome_bam        = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bai        = FASTQ_ALIGN_STAR.out.bai
        ch_transcriptome_bam = FASTQ_ALIGN_STAR.out.bam_transcript
        ch_star_log          = FASTQ_ALIGN_STAR.out.log_final
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})
    }

    if(params.with_umi) {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS (
            ch_genome_bam.join(ch_genome_bai, by: [0]),
            params.umitools_dedup_stats
        )
        ch_genome_bam = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam
        ch_genome_bai = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bai
        ch_samtools_stats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.stats
        ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.flagstat
        ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.idxstats

        ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.versions)
    }

    QUALITY_CONTROL (
        ch_genome_bam,
        PREPARE_GENOME.out.gene_bed
    )
    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)

    COVERAGE_GRAPHS (
        ch_genome_bam,
        ch_genome_bai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai
    )
    ch_versions = ch_versions.mix(COVERAGE_GRAPHS.out.versions)

    //
    // SUBWORKFLOW: Transcript indetification
    //
    ch_genome_bam.map {
        meta, bam ->
        fmeta = meta.findAll { it.key != 'read_group' }
        // Split and take the first element
        fmeta.id = fmeta.id.split('_')[0]
        [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_sort_bam }

    TRANSCRIPT_INDENTIFICATION (
        ch_sort_bam,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.chrom_sizes,
    )
    ch_grohmm_multiqc = TRANSCRIPT_INDENTIFICATION.out.grohmm_td_plot.collect()
    ch_homer_multiqc = TRANSCRIPT_INDENTIFICATION.out.homer_peaks
    ch_homer_multiqc = ch_homer_multiqc.mix(TRANSCRIPT_INDENTIFICATION.out.homer_tagdir)
    ch_versions = ch_versions.mix(TRANSCRIPT_INDENTIFICATION.out.versions)

    SUBREAD_FEATURECOUNTS_PREDICTED (
        ch_sort_bam.combine(
            BED2SAF (
                TRANSCRIPT_INDENTIFICATION.out.transcript_beds
            ).saf.map { it[1] }
        )
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_PREDICTED.out.versions.first())

    SUBREAD_FEATURECOUNTS_GENE (
        ch_sort_bam.combine(PREPARE_GENOME.out.gtf)
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_GENE.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.preseq_ccurve.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.preseq_lcextrap.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.readdistribution_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.readduplication_seq_xls.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.readduplication_pos_xls.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(QUALITY_CONTROL.out.inferexperiment_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_grohmm_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_homer_multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS_PREDICTED.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS_GENE.out.summary.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
