process GROHMM_TRANSCRIPTCALLING {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' :
    //     'quay.io/biocontainers/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' }"

    input:
    tuple val(meta), path(bams), path(bais), path(tuning)
    path gtf

    output:
    tuple val(meta), path("*.transcripts.txt"), emit: transcripts
    tuple val(meta), path("*.eval.txt")       , emit: eval
    tuple val(meta), path("*.transcripts.bed"), emit: transcripts_bed
    tuple val(meta), path("*.tdFinal.txt")    , emit: td
    tuple val(meta), path("*.tdplot_mqc.jpg") , emit: td_plot
    path  "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    transcriptcalling_grohmm.R \\
        --bam_file ${bams} \\
        --tuning_file ${tuning_file} \\
        --outprefix ${prefix} \\
        --gtf $gtf \\
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
