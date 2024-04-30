process SUBREAD2BED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas==1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2':
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(subread)
    
    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions
    
    script:
    template "subread2bed.py"
}