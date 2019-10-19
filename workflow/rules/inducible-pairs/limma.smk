rule genes_limma:
    """
    Runes the gene counts through limma script
    """
    input:
        "results/2019-06-03/{cell}/counts/merged.txt"
    output:
        # "results/2019-06-26/genes/limma/{cell}_limma.txt",
        "results/2019-06-26/genes/limma/{cell}_fig1_limma.png",
        "results/2019-06-26/genes/limma/{cell}_fig2_limma.png",
        "results/2019-06-26/genes/limma/{cell}_fig3_limma.png",
        # "results/2019-06-26/genes/limma/{cell}_fig4_limma.png", # FIXME
        # "results/2019-06-26/genes/limma/{cell}_results.txt",
    params:
    conda:
        "../../envs/edgeR.yaml"
    threads: 4
    script:
        "../../scripts/dge.R"

rule eRNA_limma:
    input:
        "results/2019-06-03/eRNA/counts/{cell}_merged.txt"
    output:
        "results/2019-06-26/eRNA/limma/{cell}_fig1_limma.png",
        "results/2019-06-26/eRNA/limma/{cell}_fig2_limma.png",
        "results/2019-06-26/eRNA/limma/{cell}_fig3_limma.png",
    params:
    conda:
        "../../envs/edgeR.yaml"
    threads: 4
    script:
        "../../scripts/dge.R"
