////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
def alignerList         = ['bwa', 'bwamem2']
if (!params.skip_alignment) {
    if (!alignerList.contains(params.aligner)) {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: ${alignerList.join(', ')}"
    }
    prepareToolIndices << params.aligner
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

def subread_featurecounts_options                      = modules['subread_featurecounts']
def subread_featurecounts_gene_options                 = subread_featurecounts_options.clone()
def subread_featurecounts_predicted_options            = subread_featurecounts_options.clone()
subread_featurecounts_predicted_options['publish_dir'] = "${params.aligner}/featurecounts/predicted"
subread_featurecounts_predicted_options.args           = " -F \"SAF\""

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']] )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
def gffread_options         = modules['gffread']
if (!params.save_reference) { gffread_options['publish_files'] = false }

def bwa_align_options            = modules['bwa_align']
// TODO
// bwa_align_options.args          += params.save_unaligned ? " --outReadsUnmapped Fastx" : ''
// if (params.save_align_intermeds)  { bwa_align_options.publish_files.put('bam','') }
// if (params.save_unaligned)        { bwa_align_options.publish_files.put('fastq.gz','unmapped') }

def samtools_sort_options = modules['samtools_sort']
// TODO
// if (['bwa'].contains(params.aligner)) {
//     if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// } else {
//     if (params.save_align_intermeds || params.skip_markduplicates) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// }

// Local: Sub-workflows
include { INPUT_CHECK           } from '../subworkflows/local/input_check'       addParams( options: [:]                          )
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'    addParams( genome_options: publish_genome_options, index_options: publish_index_options, gffread_options: gffread_options )
include { ALIGN_BWA             } from '../subworkflows/local/align_bwa'         addParams( align_options: bwa_align_options, samtools_options: samtools_sort_options )
include { ALIGN_BWAMEM2         } from '../subworkflows/local/align_bwamem2'     addParams( align_options: bwa_align_options, samtools_options: samtools_sort_options )
include { HOMER_GROSEQ          } from '../subworkflows/nf-core/homer_groseq.nf' addParams( options: [:]                          )
include { GROHMM                } from '../subworkflows/local/grohmm'            addParams( options: [:]                          )
// nf-core/modules: Modules
include { FASTQC                                                   } from '../modules/nf-core/software/fastqc/main'                addParams( options: modules['fastqc']                       )
include { CAT_FASTQ                                                } from '../modules/nf-core/software/cat/fastq/main'             addParams( options: cat_fastq_options                )
include { BED2SAF                                                  } from '../modules/local/bed2saf'                       addParams(                                                  )
include { PICARD_MERGESAMFILES                                     } from '../modules/nf-core/software/picard/mergesamfiles/main'        addParams( options: modules['picard_mergesamfiles'] )
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_PREDICTED } from '../modules/nf-core/software/subread/featurecounts/main' addParams( options: subread_featurecounts_predicted_options )
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE      } from '../modules/nf-core/software/subread/featurecounts/main' addParams( options: subread_featurecounts_options           )
include { MULTIQC                                                  } from '../modules/nf-core/software/multiqc/main'               addParams( options: multiqc_options                         )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow GROSEQ {
    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        prepareToolIndices
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
    )
    .map {
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

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    /*
     * MODULE: Run FastQC
     */
    FASTQC (
        ch_cat_fastq
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    

    /*
     * SUBWORKFLOW: Alignment with BWA
     */
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
        ch_software_versions = ch_software_versions.mix(ALIGN_BWA.out.bwa_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BWA.out.samtools_version.first().ifEmpty(null))
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
        ch_software_versions = ch_software_versions.mix(ALIGN_BWAMEM2.out.bwa_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BWAMEM2.out.samtools_version.first().ifEmpty(null))
    }

    ch_genome_bam.map {
        meta, bam ->
        fmeta = meta.findAll { it.key != 'read_group' }
        fmeta.id = "meta"
        [ fmeta, bam ] }
        .groupTuple(by: [0])
        .map { it ->  [ it[0], it[1].flatten() ] }
        .set { meta_bam }

    PICARD_MERGESAMFILES (
        meta_bam
    )

    ch_homer_multiqc = Channel.empty()
    if (params.transcript_identification == 'grohmm') {
        /*
         * SUBWORKFLOW: Transcript indetification with GROHMM
         */
        GROHMM ( PICARD_MERGESAMFILES.out.bam )

        SUBREAD_FEATURECOUNTS_PREDICTED (
            ch_genome_bam.combine( BED2SAF ( GROHMM.out.bed ) )
        )
    } else if (params.transcript_identification == 'homer') {
        /*
         * SUBWORKFLOW: Transcript indetification with homer
         */
        HOMER_GROSEQ(
            ch_genome_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_homer_multiqc = HOMER_GROSEQ.out.tag_dir
    }

    SUBREAD_FEATURECOUNTS_GENE (
        ch_genome_bam.combine(PREPARE_GENOME.out.gtf)
    )
    ch_software_versions = ch_software_versions.mix(SUBREAD_FEATURECOUNTS_GENE.out.version.first().ifEmpty(null))

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */  
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_homer_multiqc.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
