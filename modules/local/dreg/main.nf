process DREG_PEAKCALLING {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "samesense/dreg-docker"

    input:
    tuple val(meta), path(plus), path(minus)
    path model

    output:
    tuple val(meta), path("*.bedGraph.gz"), emit: bedgraph
    path  "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bash /dREG/run_peakcalling.bsh \\
        $plus \\
        $minus \\
        $prefix \\
        $model \\
        $task.cpus \\
        0
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        dREG: \$(Rscript -e "library(dREG); cat(as.character(packageVersion('dREG')))")
    END_VERSIONS
    """
}
