rule alignHistones:
    """
    Aligns the histones for a cell line using bowtie2
    TODO Use star
    TODO Use an Input function
    """
    input:
        sample=["data/2018-11-13/{cell}_{unit}.fastq"]
    output:
        "data/2018-11-13/{genome}/{cell}_{unit}.bam"
    log:
        "logs/{genome}/{cell}/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/{genome}/genome",
        extra=""
    wrapper:
        "0.35.1/bio/bowtie2/align"

rule convert_Histones_to_Bed:
    """
    Coverts Histone Bam files to BED format
    """
    input:
        "data/2018-11-13/{genome}/{cell}_{unit}.bam",
    output:
        "data/2018-11-13/{genome}/{cell}_{unit}.bed",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bamToBed -i {input} | sortBed -i - > {output}"

rule eRNAs:
    """
    Takes in units with genes removed
    Keeps anything that intersects with H3K27ac or H3K4me1
    """
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
    shell:
        "bedtools intersect -a {input.no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"

# rule eRNAs_samples_GM:
#     """
#     Takes in units with genes removed
#     Keeps anything that intersects with H3K27ac or H3K4me1
#     TODO
#     """
#     input:
#         no_genes="results/2018-11-09/{genome}/{sample}-{unit}_transcripts_noGenes.bed",
#         H3K27ac="data/2018-11-13/{genome}/GM_H3K27ac.bed",
#         H3K4me1="data/2018-11-13/{genome}/GM_H3K4me1.bed",
#     output:
#         "results/2018-12-02/{genome}/.bed"
#     log:
#         "logs/{genome}/{sample}-{unit}/HistonesIntersect.log"
#     conda:
#         "../../envs/bedtools.yaml"
#     shell:
#         "bedtools intersect -a {input.no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"

# rule eRNAs_0h_IMR:
#     """
#     Takes in units with genes removed
#     Keeps anything that intersects with H3K27ac or H3K4me1
#     TODO
#     """
#     input:
#         no_genes="results/2018-11-09/{genome}/{sample}-{unit}_transcripts_noGenes.bed",
#         H3K27ac="data/2018-11-13/{genome}/IMR_H3K27ac.bed",
#         H3K4me1="data/2018-11-13/{genome}/IMR_H3K4me1.bed",
#     output:
#         "results/2018-12-02/{genome}/IMR_eRNA_0h.bed"
#     log:
#         "logs/{genome}/IMR/HistonesIntersect.log"
#     conda:
#         "../../envs/bedtools.yaml"
#     shell:
#         "bedtools intersect -a {input.no_genes} -b {input.H3K27ac} {input.H3K4me1} -sorted -u -bed > {output} 2> {log}"
