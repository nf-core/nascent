rule IMR_genes_countReads:
    input:
        bams=expand("results/2018-12-01/IMR/{unit}.bam",unit=IMR_SAMPLES),
        genes="data/2018-11-28/genes.gtf",
        # genes="results/2019-02-28/GM19_meta_Genes.gtf",
    output:
        "results/2019-03-03/IMR_counts.txt"
    log:
        "logs/GM19/countReads.log"
    threads: 4
    conda:
        "../../envs/subread.yaml"
    shell:
        "featureCounts -T {threads} -t exon -g gene_id -a {input.genes} \
        -o {output} {input.bams}"
