process GROHMM_TRANSCRIPTCALLING {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/grohmm:03357458e0821bcb' :
        'community.wave.seqera.io/library/grohmm:833aa94cad4202ac' }"

    input:
    tuple val(meta), path(bams), path(bais), path(tuning_file)
    path gtf

    output:
    tuple val(meta), path("*.transcripts.txt"), emit: transcripts
    tuple val(meta), path("*.eval.txt")       , emit: eval
    tuple val(meta), path("*.transcripts.bed"), emit: transcripts_bed
    tuple val(meta), path("*.tdFinal.txt")    , emit: td
    tuple val(meta), path("*.tdplot_mqc.png") , emit: td_plot
    tuple val(meta), path("*.tdFinal_mqc.csv") , emit: mqc_csv
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
        --memory ${task.memory.toMega()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-grohmm: \$(Rscript -e "library(groHMM); cat(as.character(packageVersion('groHMM')))")
    END_VERSIONS
    """
}
