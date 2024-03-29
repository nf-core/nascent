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
    tuple val(meta), path("${prefix}.pos.bw"), path("${prefix}.neg.bw"), emit: dreg_bigwig
    tuple val(meta), path("${prefix}.bedGraph"), emit: dreg_bg

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        template "proseq2.0_se"
    } else {
        template "proseq2.0_pe"
    }
}
