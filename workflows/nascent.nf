/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowNascent.initialise(params, log)

// HACK Rework this because of nf-validation
def prepareToolIndices = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BED2SAF } from '../modules/local/bed2saf'

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_BWAMEM2 } from '../subworkflows/local/align_bwamem2/main'
include { ALIGN_DRAGMAP } from '../subworkflows/local/align_dragmap/main'
include { QUALITY_CONTROL } from '../subworkflows/local/quality_control.nf'
include { COVERAGE_GRAPHS } from '../subworkflows/local/coverage_graphs.nf'
include { TRANSCRIPT_INDENTIFICATION } from '../subworkflows/local/transcript_identification.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from '../modules/nf-core/fastqc/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main'
include {
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_PREDICTED } from '../modules/nf-core/subread/featurecounts/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_ALIGN_BWA } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { FASTQ_ALIGN_BOWTIE2 } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NASCENT {

    ch_versions = Channel.empty()
    ch_nascent_logo = Channel.fromPath("$projectDir/docs/images/nf-core-nascent_logo_light.png")

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        prepareToolIndices,
        params.fasta,
        params.gtf,
        params.gff,
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.first())
    ch_fasta = PREPARE_GENOME.out.fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }

    //
    // Create input channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
            } else {
                return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
            }
        }
        .groupTuple()
        .map {
            WorkflowNascent.validateInput(it)
        }
        .map {
            meta, fastqs ->
            return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ch_reads = Channel.empty()
    if(!params.skip_trimming) {
        FASTP ( ch_fastq, [], false, false )
        ch_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    } else {
        ch_reads = ch_fastq
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

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary = WorkflowNascent.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNascent.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        .mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        .mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        .mix(ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]))
        .mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
        .mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
        .mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.preseq_ccurve.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.preseq_lcextrap.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.readdistribution_txt.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.readduplication_seq_xls.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.readduplication_pos_xls.collect{it[1]}.ifEmpty([]))
        .mix(QUALITY_CONTROL.out.inferexperiment_txt.collect{it[1]}.ifEmpty([]))
        .mix(ch_grohmm_multiqc.collect{it[1]}.ifEmpty([]))
        .mix(ch_homer_multiqc.collect{it[1]}.ifEmpty([]))
        .mix(SUBREAD_FEATURECOUNTS_PREDICTED.out.summary.collect{it[1]}.ifEmpty([]))
        .mix(SUBREAD_FEATURECOUNTS_GENE.out.summary.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
