def VERSION = '4.11' // Version information not provided by tool on CLI

process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e7cf478db279450ded3321ece4ab34e14ee09fff:24fc913fe7d89b2142a4adaf04f4b934589498b4-0' :
        'quay.io/biocontainers/mulled-v2-e7cf478db279450ded3321ece4ab34e14ee09fff:24fc913fe7d89b2142a4adaf04f4b934589498b4-0' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_tagdir"), emit: tagdir
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeTagDirectory \\
        ${prefix}_tagdir \\
        $bam \\
        -genome $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
