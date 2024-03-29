process DREG_PREP {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:4238fb992d2a93e648108c86e3a9f51348e834a9-0' :
        'biocontainers/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:4238fb992d2a93e648108c86e3a9f51348e834a9-0' }"

    input:
    tuple val(meta), path(bam_file), val(index)
    path  sizes
    val assay_type

    output:
    tuple val(meta), path("${prefix}_plus.rpm.bw"), emit: plus_rpm_bigwig
    tuple val(meta), path("${prefix}_minus.rpm.bw"), emit: minus_rpm_bigwig
    tuple val(meta), path("${prefix}_plus.bw"), emit: plus_bigwig
    tuple val(meta), path("${prefix}_minus.bw"), emit: minus_bigwig

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "proseq2.0"
}
