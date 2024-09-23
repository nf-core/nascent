process GROHMM_TRANSCRIPTCALLING {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/grohmm:a660d9c3942c9b85' :
    'community.wave.seqera.io/library/grohmm:780b8693bdaa87b9' }"

    input:
    tuple val(meta), path(bams), path(bais), path(tuning_file)
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
