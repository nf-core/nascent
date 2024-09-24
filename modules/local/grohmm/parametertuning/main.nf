process GROHMM_PARAMETERTUNING {
    tag "$meta.id|$UTS|$LtProbB"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/grohmm:a660d9c3942c9b85' :
        'community.wave.seqera.io/library/grohmm:780b8693bdaa87b9' }"

    input:
    tuple val(meta), path(bams), path(bais), val(UTS), val(LtProbB)
    path gtf

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
