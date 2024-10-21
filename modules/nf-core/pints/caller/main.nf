process PINTS_CALLER {
    tag "$meta.id" + "${chr_name ? ' | ' + chr_name : ''}"
    label 'process_high'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypints:1.1.8--pyh7cba7a3_0' :
        'biocontainers/pypints:1.1.8--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    val assay_type
    each chr_name // optional

    output:
    tuple val(meta), path("*_divergent_peaks.bed")     , optional:true, emit: divergent_TREs
    tuple val(meta), path("*_bidirectional_peaks.bed") , optional:true, emit: bidirectional_TREs
    tuple val(meta), path("*_unidirectional_peaks.bed"), optional:true, emit: unidirectional_TREs
    tuple val(meta), path("peakcalling_*.log")                        , emit: peakcalling_log
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" + (chr_name ? '_' + chr_name : '_all')
    def chr = chr_name ? "--chromosome-start-with $chr_name" : ''
    // TODO handle bigwigs
    // def input_type  = ("${input[0]}".endsWith(".bam")) ? "--bam-file $input" :
    //                    ("$input".contains(".bw")) ? "--bw-pl ${input[0]} --bw-mn ${input[1]}" :
    //                    error "Please use bam or BigWig files"
    """
    pints_caller \\
        --bam-file $bams \\
        --save-to . \\
        --file-prefix $prefix \\
        --thread $task.cpus \\
        --dont-check-updates \\
        --exp-type $assay_type \\
        $chr \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_caller --version)
    END_VERSIONS
    """
}
