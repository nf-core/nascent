rule alignHistones_hg18:
    input:
        sample=["data/2018-11-13/GM_{unit}.fastq"]
    output:
        "results/2018-11-13/GM_{unit}_hg18.bam"
    log:
        "logs/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule convert_Histones_to_Bed:
    input:
        "results/2018-11-13/GM_{unit}_hg18.bam",
    output:
        "results/2018-11-24/GM_{unit}_hg18.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} > {output}"

rule sort_Histones_Bed:
    input:
        "results/2018-11-24/GM_{unit}_hg18.bed",
    output:
        "results/2018-11-24/GM_{unit}_hg18.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 4
    shell:
        "sortBed -i {input} > {output}"

rule HistonesIntersect:
    input:
        GM_no_genes="results/2018-11-09/groseq_GM_noGenes.bed",
        H3K27ac="results/2018-11-24/GM_H3K27ac_hg18.sorted.bed",
        H3K4me1="results/2018-11-24/GM_H3K4me1_hg18.sorted.bed",
    output:
        "results/2018-11-10/eRNA_GM_hg18.bed"
    log:
        "logs/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.GM_no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
