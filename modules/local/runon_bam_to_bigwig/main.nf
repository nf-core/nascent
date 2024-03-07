process RUNON_BAM_TO_BIGWIG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/pipelines/nascent/runon_bam_to_bigwig'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)
    path(chrom_info)

    output:
    tuple val(meta), path("*.bigWig")   , emit: bigwig, optional: true
    tuple val(meta), path("*.bedgraph") , emit: bedgraph, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '-SE' : '-PE'

    """
    RunOnBamToBigWig.bsh \\
        --thread=${task.cpus} \\
        $paired_end \\
        -c $chrom_info \\



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        # TODO
        # samtools (version: 1.9 http://www.htslib.org/download/)
        # bedtools v2.28.0 (http://bedtools.readthedocs.org/en/latest/)
        # bedGraphToBigWig (from the Kent source utilities http://hgdownload.cse.ucsc.edu/admin/exe/)
    END_VERSIONS
    """
}
