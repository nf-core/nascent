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
    // template "proseq2.0"
    """
    bamToBed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
        awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
        > ${prefix}.dreg.bed
    sortBed -i ${prefix}.dreg.bed > ${prefix}.dreg.sort.bed

    echo "positive strand processed to bedGraph"

    bedtools genomecov \
        -bg \
        -i ${prefix}.dreg.sort.bed \
        -g ${params.chrom_sizes} \
        -strand + \
        > ${prefix}.pos.bedGraph

    sortBed \
        -i ${prefix}.pos.bedGraph \
        > ${prefix}.pos.sort.bedGraph

    bedtools genomecov \
        -bg \
        -i ${prefix}.dreg.sort.bed \
        -g ${params.chrom_sizes} \
        -strand - \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${prefix}.neg.bedGraph

    sortBed \
        -i ${prefix}.neg.bedGraph \
        > ${prefix}.neg.sort.bedGraph

    echo "negative strand processed to bedGraph"

    ${params.bedGraphToBigWig} ${prefix}.pos.sort.bedGraph ${params.chrom_sizes} ${prefix}.pos.bw
    ${params.bedGraphToBigWig} ${prefix}.neg.sort.bedGraph ${params.chrom_sizes} ${prefix}.neg.bw

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph
    """
}
