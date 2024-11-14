process HOMER_CREATEUNIQMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(mappable_regions)

    output:
    tuple val(meta), path("uniqmap/"), emit: uniqmap_dir
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11'
    """
    mkdir -p uniqmap
    homerTools special uniqmap uniqmap/ $mappable_regions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '4.11'
    """
    mkdir -p uniqmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
