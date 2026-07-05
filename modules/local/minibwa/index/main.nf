process MINIBWA_INDEX {
    tag "$fasta"
    label 'process_high'
    memory { 18.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minibwa_htslib_samtools:8111e1ae64c560b2' :
        'community.wave.seqera.io/library/minibwa_htslib_samtools:6d867025b2559348' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("minibwa"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args = task.ext.args ?: ''
    """
    mkdir minibwa
    minibwa \\
        index \\
        $args \\
        -t $task.cpus \\
        $fasta \\
        minibwa/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir minibwa
    touch minibwa/${prefix}.l2b
    touch minibwa/${prefix}.mbw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
