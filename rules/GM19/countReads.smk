rule GM19_genes_featureCounts:
    input:
        bams=expand("results/2018-11-27/GM19/{unit}.bam",unit=GM_SAMPLES),
    output:
        "results/2019-02-28/GM19_gene_counts.rds"
    log:
        "logs/GM19/genes/countReads.log"
    threads: 4
    conda:
        "../../envs/subread.yaml"
    script:
        "../../scripts/featureCount.R"

rule GM19_genes_count_cleanup:
    input:
        "results/2019-02-28/GM19_gene_counts.txt"
    output:
        "results/2019-02-28/GM19_gene_counts_clean.txt"
    shell:
        "tail -n +1 {input} > {output}"

rule GM19_eRNAs_countReads:
    input:
        bams=expand("results/2018-11-27/GM19/{unit}.bam",unit=GM_SAMPLES),
        eRNAs="results/2019-01-31/GM19_eRNA_peak.gff",
    output:
        "results/2019-02-28/GM19_eRNA_counts.txt"
    log:
        "logs/GM19/eRNA_countReads.log"
    threads: 4
    conda:
        "../../envs/subread.yaml"
    shell:
        "featureCounts -T {threads} -a {input.eRNAs}\
        -F GFF -o {output} {input.bams}"
