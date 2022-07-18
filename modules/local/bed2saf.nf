process BED2SAF {
    tag "$bed"
    label 'process_low'

    conda (params.enable_conda ? "pandas" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' :
        'quay.io/biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    path bed

    output:
    path '*.saf' , emit: saf

    script: // This script is bundled with the pipeline, in nf-core/nascent/bin/
    """
    bed2saf.py -o ${bed.baseName}.saf $bed
    """
}
