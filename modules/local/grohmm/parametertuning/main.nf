process GROHMM_PARAMETERTUNING {
    tag "$meta.id|$UTS|$LtProbB"
    label 'process_high'
    // array 10

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b929af5662486ba6ce2d27eb501e5c7ec71ca7dd8e333fe5d3dcf2803d87cf67/data' :
        'community.wave.seqera.io/library/grohmm:833aa94cad4202ac' }"

    input:
    tuple val(meta), path(bams), path(bais), val(UTS), val(LtProbB)
    path gxf

    output:
    tuple val(meta), path("*.tuning.csv"), emit: tuning
    tuple val(meta), path("*.tuning.consensus.bed"), emit: bed
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${UTS}_${LtProbB}"
    """
    grohmm_parametertuning.R \\
        --bam_file ${bams} \\
        --outprefix ${prefix} \\
        --gxf $gxf \\
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
