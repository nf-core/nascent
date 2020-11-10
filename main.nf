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

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

include { CAT_FASTQ } from './modules/local/process/cat_fastq'                   addParams( options: cat_fastq_options                                               ) 

include { INPUT_CHECK     } from './modules/local/subworkflow/input_check'     addParams( options: [:] )
include { PREPARE_GENOME  } from './modules/local/subworkflow/prepare_genome'  addParams( gffread_options: gffread_options, genome_options: publish_genome_options )


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
}
