#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/groseq
========================================================================================
 nf-core/groseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/groseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run nf-core/groseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}


////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

include { CAT_FASTQ } from './modules/local/process/cat_fastq'                   addParams( options: cat_fastq_options                                               ) 


/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
def gffread_options         = modules['gffread']
if (!params.save_reference) { gffread_options['publish_files'] = false }

def bwa_index_options = modules['bwa_index']
if (!params.save_reference)     { bwa_index_options['publish_files'] = false }

def bwa_mem_options            = modules['bwa_mem']
if (params.save_align_intermeds)  { bwa_mem_options.publish_files.put('bam','') }
if (params.save_unaligned)        { bwa_mem_options.publish_files.put('fastq.gz','unmapped') }

def samtools_sort_options = modules['samtools_sort']
if (['star','hisat2'].contains(params.aligner)) {
    if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
} else {
    if (params.save_align_intermeds || params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}

include { INPUT_CHECK    } from './modules/local/subworkflow/input_check'    addParams( options: [:] )
include { PREPARE_GENOME } from './modules/local/subworkflow/prepare_genome' addParams( gffread_options: gffread_options, genome_options: publish_genome_options )
include { ALIGN_BWA      } from './modules/local/subworkflow/align_bwa'      addParams( index_options: bwa_index_options, align_options: bwa_mem_options, samtools_options: samtools_sort_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
def umitools_extract_options    = modules['umitools_extract']
umitools_extract_options.args  += params.umitools_extract_method ? " --extract-method=${params.umitools_extract_method}" : ''
umitools_extract_options.args  += params.umitools_bc_pattern     ? " --bc-pattern='${params.umitools_bc_pattern}'"       : ''
if (params.save_umi_intermeds)  { umitools_extract_options.publish_files.put('fastq.gz','') }

def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

include { FASTQC_UMITOOLS_TRIMGALORE } from './modules/nf-core/subworkflow/fastqc_umitools_trimgalore' addParams( fastqc_options: modules['fastqc'], umitools_options: umitools_extract_options, trimgalore_options: trimgalore_options               )

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


workflow {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.gene_bed,
        params.additional_fasta
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, bam ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, bam ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_cat_fastq }

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ (
        ch_cat_fastq
    )

    /*
     * SUBWORKFLOW: Read QC, extract UMI and trim adapters
     */
    FASTQC_UMITOOLS_TRIMGALORE (
        CAT_FASTQ.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.umitools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    ch_trimmed_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads

    /*
     * SUBWORKFLOW: Alignment with bwa
     */
    ch_genome_bam        = Channel.empty()
    ch_genome_bai        = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_star_multiqc      = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bwa') {
        ALIGN_BWA (
            ch_trimmed_reads,
            params.star_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_genome_bam        = ALIGN_BWA.out.bam
        ch_genome_bai        = ALIGN_BWA.out.bai
        ch_samtools_stats    = ALIGN_BWA.out.stats
        ch_samtools_flagstat = ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA.out.idxstats
        ch_software_versions = ch_software_versions.mix(ALIGN_BWA.out.bwa_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BWA.out.samtools_version.first().ifEmpty(null))
    }
}
