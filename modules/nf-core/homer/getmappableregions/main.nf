process HOMER_GETMAPPABLEREGIONS {
    tag "${fasta_files[0].baseName}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'community.wave.seqera.io/library/homer:4.11--e7dc4a041f589d54' }"

    input:
    path(fasta_files, arity: '1..*')
    val(parallel_sequences)
    val(read_length)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta_files[0].baseName}"
    def VERSION = '4.11'
    """
    getMappableRegions \\
        $parallel_sequences \\
        $read_length \\
        $fasta_files \\
        > ${prefix}.${read_length}nt.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11'
    """
    touch ${prefix}.${read_length}nt.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
