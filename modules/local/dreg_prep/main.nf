process DREG_PREP {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nf-core/nascent/dreg_prep:dreg_prep--18e6e4570a433d65' :
        'nf-core/nascent/dreg_prep:b6f1c5d2a964415f' }"

    input:
    tuple val(meta), path(bam_file), val(index)
    path  sizes

    output:
    tuple val(meta), path("${prefix}.pos.bw"), path("${prefix}.neg.bw"), emit: dreg_bigwig
    tuple val(meta), path("${prefix}.bedGraph"), emit: dreg_bg

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo "Creating BigWigs suitable as inputs to dREG"

        bamToBed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
            awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
            > ${prefix}.dreg.bed
        sortBed -i ${prefix}.dreg.bed > ${prefix}.dreg.sort.bed

        echo "positive strand processed to bedGraph"

        bedtools genomecov \
            -bg \
            -i ${prefix}.dreg.sort.bed \
            -g ${sizes} \
            -strand + \
            > ${prefix}.pos.bedGraph

        sortBed \
            -i ${prefix}.pos.bedGraph \
            > ${prefix}.pos.sort.bedGraph

        bedtools genomecov \
            -bg \
            -i ${prefix}.dreg.sort.bed \
            -g ${sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${prefix}.neg.bedGraph

        sortBed \
            -i ${prefix}.neg.bedGraph \
            > ${prefix}.neg.sort.bedGraph

        echo "negative strand processed to bedGraph"

        bedGraphToBigWig ${prefix}.pos.sort.bedGraph ${sizes} ${prefix}.pos.bw
        bedGraphToBigWig ${prefix}.neg.sort.bedGraph ${sizes} ${prefix}.neg.bw

        cat ${prefix}.pos.bedGraph \
            ${prefix}.neg.bedGraph \
            > ${prefix}.unsorted.bedGraph

        sortBed \
            -i ${prefix}.unsorted.bedGraph \
            > ${prefix}.bedGraph

        echo "bedGraph to bigwig done"
        """
    } else {
        if (params.forwardStranded) {
            """
            samtools view -@ 16 -bf 0x2 ${bam_file} | samtools sort -n -@ 16 \
                > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
                | awk 'BEGIN{OFS="\t"} (\$9 == "+") {print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9 == "-") {print \$1,\$3-1,\$3,\$7,\$8,\$9}' \
                | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${sizes} \
                -strand + \
                > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${sizes} \
                -strand - \
                > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
                | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
                > ${prefix}.neg.bedGraph

            bedGraphToBigWig ${prefix}.pos.bedGraph \
                ${sizes} ${prefix}.pos.bw

            bedGraphToBigWig ${prefix}.neg.bedGraph \
                ${sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
                ${prefix}.neg.bedGraph \
                > ${prefix}.unsorted.bedGraph

            sortBed \
                -i ${prefix}.unsorted.bedGraph \
                > ${prefix}.bedGraph
            """
        } else {
            """
            samtools view -@ 16 -bf 0x2 ${bam_file} | samtools sort -n -@ 16 \
                > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
                | awk 'BEGIN{OFS="\t"} (\$10 == "+") {print \$1,\$5,\$5+1,\$7,\$8,\$10}; (\$10 == "-") {print \$1,\$6-1,\$6,\$7,\$8,\$10}' \
                | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${sizes} \
                -strand + \
                > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${sizes} \
                -strand - \
                > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
                | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
                > ${prefix}.neg.bedGraph

            bedGraphToBigWig ${prefix}.pos.bedGraph \
                ${sizes} ${prefix}.pos.bw

            bedGraphToBigWig ${prefix}.neg.bedGraph \
                ${sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
                ${prefix}.neg.bedGraph \
                > ${prefix}.unsorted.bedGraph

            sortBed \
                -i ${prefix}.unsorted.bedGraph \
                > ${prefix}.bedGraph
            """
        }
    }
}
