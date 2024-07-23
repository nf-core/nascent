process PINTS_CALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypints:1.1.8--pyh7cba7a3_0' :
        'biocontainers/pypints:1.1.8--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(bams), path(bais)
    val assay_type

    output:
    tuple val(meta), path("*_divergent_peaks.bed")     , optional:true, emit: divergent_TREs
    tuple val(meta), path("*_bidirectional_peaks.bed") , optional:true, emit: bidirectional_TREs
    tuple val(meta), path("*_unidirectional_peaks.bed"), optional:true, emit: unidirectional_TREs
    tuple val(meta), path("peakcalling_*.log")                        , emit: peakcalling_log

    eval("pints_caller --version")               , topic: version
    eval("python --version | sed 's/Python //g'"), topic: version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
        $args
   """
}
