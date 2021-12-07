/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
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
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK           } from '../subworkflows/local/input_check'
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'
include { GROHMM                } from '../subworkflows/local/grohmm'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { FASTQC                                                  } from '../modules/nf-core/modules/fastqc/main'
include { CAT_FASTQ                                               } from '../modules/nf-core/modules/cat/fastq/main'
include { BED2SAF                                                 } from '../modules/local/bed2saf'
include { PICARD_MERGESAMFILES                                    } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE
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
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
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

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions.first())

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

    ch_genome_bam.map {
        meta, bam ->
        fmeta = meta.findAll { it.key != 'read_group' }
        fmeta.id = fmeta.group
        [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { group_bam }

    PICARD_MERGESAMFILES (
        group_bam
    )

    ch_homer_multiqc = Channel.empty()
    if (params.transcript_identification == 'grohmm') {
        //
        // SUBWORKFLOW: Transcript indetification with GROHMM
        //
        GROHMM (
            PICARD_MERGESAMFILES.out.bam,
            PREPARE_GENOME.out.gtf
        )

        SUBREAD_FEATURECOUNTS_PREDICTED (
            ch_genome_bam.combine( BED2SAF ( GROHMM.out.bed ) )
        )
    } else if (params.transcript_identification == 'homer') {
        /*
        * SUBWORKFLOW: Transcript indetification with homer
        */
        HOMER_GROSEQ (
            PICARD_MERGESAMFILES.out.bam,
            PREPARE_GENOME.out.fasta
        )
        ch_homer_multiqc = HOMER_GROSEQ.out.tag_dir
    }

    SUBREAD_FEATURECOUNTS_GENE (
        ch_genome_bam.combine(PREPARE_GENOME.out.gtf)
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_GENE.out.versions.first())

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNascent.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_homer_multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_versions = ch_versions.mix(MULTIQC.out.versions.first())


    //
    // MODULE: Dump software versions for all tools used in the workflow
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/