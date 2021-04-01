// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Convert BED file to SAF format
 */
process BED2SAF {
    tag "$bed"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', publish_id:'') }

    conda (params.enable_conda ? "pandas" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    }

    input:
    path bed
    
    output:
    path '*.saf' , emit: saf

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    bed2saf.py -o ${bed.baseName}.saf $bed
    """
}
