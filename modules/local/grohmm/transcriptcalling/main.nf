process GROHMM_TRANSCRIPTCALLING{
    tag "$meta.id"
    label 'process_high'

    conda    (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-optparse=1.6.6 conda-forge::r-rmariadb=1.1.0 bioconda::bioconductor-edger=3.32.0 bioconda::bioconductor-grohmm=1.24.0 bioconda::bioconductor-org.hs.eg.db=3.12.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2"  : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-598ad5976af555e8404016e81a380911377ebb95:79c659213212873c7f35071e27b6732b6de63a77-0' :
        'quay.io/biocontainers/mulled-v2-598ad5976af555e8404016e81a380911377ebb95:79c659213212873c7f35071e27b6732b6de63a77-0' }"

    input:
    tuple val(meta), path(bam)
    path gtf
    path tuning

    output:
    path "*.transcripts.txt" , emit: transcripts
    path "*.eval.txt"        , emit: eval
    path "*.transcripts.bed" , emit: transcripts_bed
    path "*.tdFinal.txt"     , emit: td
    path "*.tdplot.jpg"      , emit: td_plot
    // FIXME path "*.RData"  , emit: rdata
    path  "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (params.with_tuning) {
        """
        transcriptcalling_grohmm.R \\
            --bam_file ${bam} \\
            --tuning_file ${tuning} \\
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
    } else {
        """
        transcriptcalling_grohmm.R \\
            --bam_file ${bam} \\
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
}
