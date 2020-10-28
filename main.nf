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
log.info "Hello"
