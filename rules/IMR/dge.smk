rule IMR_genes_featureCounts:
    input:
        bams=expand("results/2018-12-01/IMR/{unit}.bam",unit=IMR_SAMPLES),
    output:
        "results/2019-02-28/IMR_gene_counts.rds"
    log:
        "logs/IMR/genes/countReads.log"
    threads: 4
    conda:
        "../../envs/subread.yaml"
    script:
        "../../scripts/featureCountIMR.R"

# TODO only spits out counts
rule IMR_genes_edgeR:
    input:
        "results/2019-02-28/IMR_gene_counts.rds"
    output:
        "results/2019-03-25/IMR_rawdata_table.csv",
        "results/2019-03-25/IMR_genelen.csv",
        "results/2019-03-25/IMR_normalized.csv",
    conda:
        "../../envs/edgeR.yaml"
    threads: 4
    script:
        "../../scripts/dge.R"
