rule alignHistones:
    input:
        sample=["data/2018-11-13/{cell}_{unit}.fastq"]
    output:
        "data/2018-11-13/{genome}/{cell}_{unit}.bam"
    log:
        "logs/{genome}/{cell}/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/{genome}/genome",
        extra=""
    threads: 4
    wrapper:
        "0.35.1/bio/bowtie2/align"

rule convert_Histones_to_Bed:
    input:
        "data/2018-11-13/{genome}/{cell}_{unit}.bam",
    output:
        "data/2018-11-13/{genome}/{cell}_{unit}.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bamToBed -i {input} | sortBed -i - > {output}"

rule histonesIntersect:
    input:
        no_genes="results/2018-11-09/{genome}/{cell}_meta_transcripts_noGenes.bed",
        H3K27ac="data/2018-11-13/{genome}/{cell}_H3K27ac.bed",
        H3K4me1="data/2018-11-13/{genome}/{cell}_H3K4me1.bed",
    output:
        "results/2018-12-02/{genome}/{cell}_eRNA.bed"
    log:
        "logs/{genome}/{cell}/HistonesIntersect.log"
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "bedtools intersect -a {input.no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
