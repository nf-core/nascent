/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['bwa', 'bwamem2']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNascent.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK           } from '../subworkflows/local/input_check'
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'
include { QUALITY_CONTROL       } from '../subworkflows/local/quality_control.nf'
include { COVERAGE_GRAPHS       } from '../subworkflows/local/coverage_graphs.nf'
include { GROHMM_MAKEUCSCFILE   } from '../modules/local/grohmm/makeucscfile/main.nf'
include { GROHMM                } from '../subworkflows/local/grohmm'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                                  } from '../modules/nf-core/modules/fastqc/main'
include { CAT_FASTQ                                               } from '../modules/nf-core/modules/cat/fastq/main'
include { DREG_PEAKCALLING                                        } from '../modules/local/dreg/main'
include { BED2SAF                                                 } from '../modules/local/bed2saf'
include {
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE
    SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_PREDICTED } from '../modules/nf-core/modules/subread/featurecounts/main'
include { MULTIQC                                                 } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                             } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'


//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { ALIGN_BWA             } from '../subworkflows/nf-core/align_bwa/main'
include { ALIGN_BWAMEM2         } from '../subworkflows/nf-core/align_bwamem2/main'
include { HOMER_GROSEQ          } from '../subworkflows/nf-core/homer_groseq.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NASCENT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.first())

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    INPUT_CHECK.out.reads.map {
        meta, fastq ->
        meta.id = meta.id.split('_')[0..-2].join('_')
        [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
            single  : fastq.size() == 1
            return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )

    CAT_FASTQ.out.reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUBWORKFLOW: Alignment with BWA
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bai                 = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bwa') {
        ALIGN_BWA(
            ch_cat_fastq,
            PREPARE_GENOME.out.bwa_index,
        )
        ch_genome_bam        = ALIGN_BWA.out.bam
        ch_genome_bai        = ALIGN_BWA.out.bai
        ch_samtools_stats    = ALIGN_BWA.out.stats
        ch_samtools_flagstat = ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA.out.idxstats

        ch_versions = ch_versions.mix(ALIGN_BWA.out.versions.first())
    } else if (!params.skip_alignment && params.aligner == 'bwamem2') {
        ALIGN_BWAMEM2(
            ch_cat_fastq,
            PREPARE_GENOME.out.bwa_index,
        )
        ch_genome_bam        = ALIGN_BWAMEM2.out.bam
        ch_genome_bai        = ALIGN_BWAMEM2.out.bai
        ch_samtools_stats    = ALIGN_BWAMEM2.out.stats
        ch_samtools_flagstat = ALIGN_BWAMEM2.out.flagstat
        ch_samtools_idxstats = ALIGN_BWAMEM2.out.idxstats

        ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions.first())
    }

    QUALITY_CONTROL (
        ch_genome_bam,
        PREPARE_GENOME.out.gene_bed
    )
    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions.first())

    COVERAGE_GRAPHS (
        ch_genome_bam,
        ch_genome_bai,
        PREPARE_GENOME.out.chrom_sizes
    )
    ch_versions = ch_versions.mix(COVERAGE_GRAPHS.out.versions.first())

    //
    // SUBWORKFLOW: Transcript indetification
    //
    ch_genome_bam.map {
        meta, bam ->
        fmeta = meta.findAll { it.key != 'read_group' }
        // FIXME Group splitting is broken
        // fmeta.id = fmeta.id.split('_')[0..-2].join('_')
        [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { ch_sort_bam }

    ch_homer_multiqc = Channel.empty()
    ch_identification_bed = Channel.empty()
    if (params.transcript_identification == 'grohmm') {
        // FIXME
        log.warn "grohmm works with small test data set but not on full files."

        GROHMM (
            ch_sort_bam,
            PREPARE_GENOME.out.gtf
        )

        ch_identification_bed = GROHMM.out.bed

    } else if (params.transcript_identification == 'dreg') {
        // FIXME
        println "Docker image is broken"

        ch_dreg_model = file("ftp://cbsuftp.tc.cornell.edu/danko/hub/dreg.models/asvm.gdm.6.6M.20170828.rdata")

        DREG_PEAKCALLING (
            COVERAGE_GRAPHS.out.plus_minus,
            ch_dreg_model
        )

        ch_identification_bed = GROHMM.out.bed
    } else if (params.transcript_identification == 'homer') {
        /*
        * SUBWORKFLOW: Transcript indetification with homer
        */
        HOMER_GROSEQ (
            ch_sort_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_versions = ch_versions.mix(HOMER_GROSEQ.out.versions.first())

        ch_homer_multiqc = HOMER_GROSEQ.out.peaks
        ch_homer_multiqc = ch_homer_multiqc.mix(HOMER_GROSEQ.out.tagdir)
        ch_identification_bed = HOMER_GROSEQ.out.bed
    }

    ch_identification_bed.map {
        meta, bed ->
        bed
    }.set { ch_transcript_bed }

    SUBREAD_FEATURECOUNTS_PREDICTED (
        ch_genome_bam.combine( BED2SAF ( ch_transcript_bed ) )
    )


    SUBREAD_FEATURECOUNTS_GENE (
        ch_genome_bam.combine(PREPARE_GENOME.out.gtf)
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_GENE.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNascent.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.preseq_ccurve.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.preseq_lcextrap.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.readdistribution_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.readduplication_seq_xls.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.readduplication_pos_xls.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.inferexperiment_txt.collect{it[1]}.ifEmpty([]))
    // FIXME ch_multiqc_files = ch_multiqc_files.mix(GROHMM.out.td_plot.collect().ifEmpty([]))
    // FIXME ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS_PREDICTED.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_homer_multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_versions = ch_versions.mix(MULTIQC.out.versions.first())


    //
    // MODULE: Dump software versions for all tools used in the workflow
    //

    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
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
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
