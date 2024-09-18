process GROHMM_PARAMETERTUNING {
    tag "$meta.id|$UTS|$LtProbB"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' :
    //     'quay.io/biocontainers/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    path gtf
    val UTS
    val LtProbB

    output:
    tuple val(meta), path("*.tuning.csv"), emit: tuning
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${UTS}_${LtProbB}"
    """
    parameter_tuning.R \\
        --bam_file ${bams} \\
        --outprefix ${prefix} \\
        --gtf $gtf \\
        --uts $UTS \\
        --ltprobb $LtProbB \\
        --outdir ./ \\
        --cores $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-grohmm: \$(Rscript -e "library(groHMM); cat(as.character(packageVersion('groHMM')))")
    END_VERSIONS
    """
}
