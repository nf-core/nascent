process DEEPTOOLS_BAMCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamCoverage \\
        -b $bam \\
        -o ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotHeatmap --version | sed -e "s/plotHeatmap //g")
    END_VERSIONS
    """
}
