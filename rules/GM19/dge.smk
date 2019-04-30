# TODO only spits out counts
rule GM19_genes_edgeR:
    input:
        "results/2019-02-28/GM19_gene_counts.rds"
    output:
        "results/2019-03-25/GM19_rawdata_table.csv",
        "results/2019-03-25/GM19_genelen.csv",
        "results/2019-03-25/GM19_normalized.txt"
    conda:
        "../../envs/edgeR.yaml"
    threads: 4
    script:
        "../../scripts/dge.R"
