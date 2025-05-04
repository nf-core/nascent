process FLIP_STRAND {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(bigwig)

    output:
    tuple val(meta), path("*.flipped.bigWig"), emit: flipped_bigwig
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Process the bigWig file to flip the strand
    cat ${bigwig} | awk '{
        if (\$4 == "+") {
            \$4 = "-"
        } else if (\$4 == "-") {
            \$4 = "+"
        }
        print \$0
    }' > ${prefix}.flipped.bigWig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
