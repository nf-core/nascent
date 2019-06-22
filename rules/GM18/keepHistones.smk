rule GM18_alignHistones:
    input:
        sample=["data/2018-11-13/GM_{unit}.fastq"]
    output:
        "results/2018-11-13/GM18_{unit}.bam"
    log:
        "logs/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule GM18_convert_Histones_to_Bed:
    input:
        "results/2018-11-13/GM18{unit}.bam",
    output:
        "results/2018-11-24/GM18{unit}.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} > {output}"

rule GM18_sort_Histones_Bed:
    input:
        "results/2018-11-24/GM18_{unit}.bed",
    output:
        "results/2018-11-24/GM18_{unit}.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 4
    shell:
        "sortBed -i {input} > {output}"

rule GM18_HistonesIntersect:
    input:
        GM_no_genes="results/2018-11-09/GM18_meta_groseq_noGenes.bed",
        H3K27ac="results/2018-11-24/GM18_H3K27ac.sorted.bed",
        H3K4me1="results/2018-11-24/GM18_H3K4me1.sorted.bed",
    output:
        "results/2018-11-10/GM18_eRNA.bed"
    log:
        "logs/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.GM_no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
