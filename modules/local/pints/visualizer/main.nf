process PINTS_VISUALIZER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // NOTE Stopped publishing at 1.1.9 https://quay.io/repository/biocontainers/pypints?tab=tags
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a9e30012e1b41baf9acd1ff94e01161138d8aa17f4e97aa32f2dc4effafcd1/data'
        : 'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6'}"

    input:
    tuple val(meta), path(bam)
    val assay_type

    output:
    tuple val(meta), path("*_pl.bw"), emit: plus_bw
    tuple val(meta), path("*_mn.bw"), emit: minus_bw
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def chr_arg = chr_name ? "--chrom ${chr_name}" : ''
    // def norm_arg = norm_factor != null ? "--norm-fact ${norm_factor}" : ''
    // def rpm_arg = rpm_normalize ? "--rpm" : ''
    // def rc_arg = reverse_complement ? "--reverse-complement" : ''
    """
    pints_visualizer \\
        --bam ${bam} \\
        --exp-type ${assay_type} \\
        --output-prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_visualizer --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${chr_name}_plus.bigwig
    touch ${prefix}_${chr_name}_minus.bigwig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pints: \$(pints_visualizer --version)
    END_VERSIONS
    """
}
