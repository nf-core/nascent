rule hg19_alignHistones:
    input:
        sample=["data/2018-11-13/GM_{unit}.fastq"]
    output:
        "data/2018-11-13/GM19/GM19_{unit}.bam"
    log:
        "logs/GM19/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/hg19/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule hg19_convert_Histones_to_Bed:
    input:
        "data/2018-11-13/GM19/GM19_{unit}.bam",
    output:
        "data/2018-11-13/GM19/GM19_{unit}.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} > {output}"

rule hg19_sort_Histones_Bed:
    input:
        "data/2018-11-13/GM19/GM19_{unit}.bed",
    output:
        "data/2018-11-13/GM19/GM19_{unit}.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "sortBed -i {input} > {output}"

rule GM_hg19_HistonesIntersect:
    input:
        GM_no_genes="results/2018-11-09/GM19_meta_groseq_noGenes.bed",
        H3K27ac="data/2018-11-13/GM19/GM19_H3K27ac.sorted.bed",
        H3K4me1="data/2018-11-13/GM19/GM19_H3K4me1.sorted.bed",
    output:
        "results/2018-12-02/GM19_eRNA.bed"
    log:
        "logs/GM19/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.GM_no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
