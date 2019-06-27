# TODO only spits out counts
rule GM19_genes_edgeR:
    input:
        "results/2019-06-03/eRNA/counts/{sample}_merged.txt"
    output:
        "results/2019-06-26/dge/limma/{sample}_limma.txt"
    conda:
        "../../envs/edgeR.yaml"
    threads: 4
    script:
        "../../scripts/dge.R"

rule GM19_genes_NOIseq:
    input:
        "results/2019-06-03/eRNA/counts/GM19_merged.txt"
    output:
        directory("results/2019-06-13/DEGS/GM19"),
    # singularity:
    #     "docker://emiller88/noiseq:0.0.2"
    conda:
        "../../envs/NOISeq.yaml"
    threads: 4
    script:
        "../../scripts/NOISeq_biomart.R"
