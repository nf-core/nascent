process MINIBWA_MAP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minibwa_htslib_samtools:8111e1ae64c560b2' :
        'community.wave.seqera.io/library/minibwa_htslib_samtools:6d867025b2559348' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val sort_bam

    output:
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.log") , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    """
    INDEX=`find -L ./ -name "*.mbw" | sed 's/\\.mbw\$//' | head -n 1`

    minibwa \\
        map \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        2> >(tee ${prefix}.minibwa.log >&2) \\
        | samtools $samtools_command $args2 ${reference} --threads $task.cpus -o ${prefix}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version 2>&1 | head -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.minibwa.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version 2>&1 | head -n 1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
