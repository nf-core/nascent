rule IMR_hg19_alignHistones:
    input:
        sample=["data/2018-12-02/IMR_{unit}.fastq"]
    output:
        "data/2018-12-02/IMR_{unit}_hg19.bam"
    log:
        "logs/bowtie2/{unit}.log"
    params:
        index="data/2018-11-27/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule IMR_hg19_convert_Histones_to_Bed:
    input:
        "data/2018-12-02/IMR_{unit}_hg19.bam",
    output:
        "data/2018-12-02/IMR_{unit}_hg19.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} > {output}"

rule IMR_hg19_sort_Histones_Bed:
    input:
        "data/2018-12-02/IMR_{unit}_hg19.bed",
    output:
        "data/2018-12-02/IMR_{unit}_hg19.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "sortBed -i {input} > {output}"

rule IMR_hg19_HistonesIntersect:
    input:
        IMR_no_genes="results/2018-12-02/IMR_hg19_groseq_noGenes.bed",
        H3K27ac="data/2018-12-02/IMR_H3K27ac_hg19.sorted.bed",
        H3K4me1="data/2018-12-02/IMR_H3K4me1_hg19.sorted.bed",
    output:
        "results/2018-12-02/eRNA_IMR_hg19.bed"
    log:
        "logs/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.IMR_no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
