process dreg_prep {
    println "[Log 5]: Generating bigwig files for dREG"

    errorStrategy 'ignore'
    tag "$prefix"
    memory '60 GB'
    cpus 16
    queue 'short'

    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (params.savebg && (filename == "${prefix}.bedGraph"))       "bedgraphs/${prefix}.fivep.bedGraph"
             else if (params.savebw && (filename.indexOf("bw") > 0))        "bigwigs/$filename"
             else null
            }

    when:
    params.dreg || params.dreg_results

    input:
    set val(prefix), file(bam_file), file(index) from bam_for_dreg

    output:
    tuple val(prefix), file("${prefix}.pos.bw"), file("${prefix}.neg.bw") into dreg_bigwig
    tuple val(prefix), file("${prefix}.bedGraph") into dreg_bg

    script:
    if (params.singleEnd) {
        """
        echo "Creating BigWigs suitable as inputs to dREG"

        export CRAM_REFERENCE=${params.genome}

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
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

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
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

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
