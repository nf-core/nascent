rule IMR_alignHistones:
    input:
        sample=["data/2018-11-13/IMR_{unit}.fastq"]
    output:
        "data/2018-11-13/IMR/IMR_{unit}.bam"
    log:
        "logs/IMR/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/hg19/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule IMR_convert_Histones_to_Bed:
    input:
        "data/2018-11-13/IMR/IMR_{unit}.bam",
    output:
        "data/2018-11-13/IMR/IMR_{unit}.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} > {output}"

rule IMR_sort_Histones_Bed:
    input:
        "data/2018-11-13/IMR/IMR_{unit}.bed",
    output:
        "data/2018-11-13/IMR/IMR_{unit}.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "sortBed -i {input} > {output}"

rule IMR_HistonesIntersect:
    input:
        IMR_no_genes="results/2018-11-09/IMR_groseq_noGenes.bed",
        H3K27ac="data/2018-11-13/IMR/IMR_H3K27ac.sorted.bed",
        H3K4me1="data/2018-11-13/IMR/IMR_H3K4me1.sorted.bed",
    output:
        "results/2018-12-02/IMR_eRNA.bed"
    log:
        "logs/IMR/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.IMR_no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
