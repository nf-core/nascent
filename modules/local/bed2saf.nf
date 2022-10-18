process BED2SAF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.saf"), emit: saf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'OFS="\\t" {print \$1"."\$2"."\$3, \$1, \$2, \$3, "."}' \\
        $bed \\
        > ${bed.baseName}.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
