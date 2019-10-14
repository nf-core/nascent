rule genes_limma:
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
        "../envs/edgeR.yaml"
    threads: 4
    script:
        "../scripts/dge.R"

rule eRNA_limma:
    input:
        "results/2019-06-03/eRNA/counts/{cell}_merged.txt"
    output:
        "results/2019-06-26/eRNA/limma/{cell}_fig1_limma.png",
        "results/2019-06-26/eRNA/limma/{cell}_fig2_limma.png",
        "results/2019-06-26/eRNA/limma/{cell}_fig3_limma.png",
    params:
    conda:
        "../envs/edgeR.yaml"
    threads: 4
    script:
        "../scripts/dge.R"

rule genes_tpm:
    input:
        "results/2019-06-03/hg19/{cell}/counts/merged.txt"
    output:
        "results/2019-06-26/dge/tpm/{cell}_tpm.txt",
    params:
    conda:
        "../envs/matplotlib.yaml"
    threads: 4
    script:
        "../scripts/tpm.py"

rule genes_foldchange:
    input:
        "results/2019-06-26/dge/tpm/{cell}_tpm.txt",
    output:
        "results/2019-06-26/dge/foldchange/{cell}_foldchange.tsv",
    params:
        # cutoff=0.5
    conda:
        "../envs/matplotlib.yaml"
    threads: 4
    script:
        "../scripts/foldchange.py"

rule genes_de_geneid:
    input:
        genes="results/2019-06-26/dge/foldchange/{cell}_foldchange.tsv",
        refseq="data/2018-11-09/hg19/genes.gtf",
    output:
        "results/2019-06-26/dge/rg/{cell}_de_geneid.txt",
    conda:
        "../envs/gawk.yaml"
    threads: 1
    shell:
        """awk -F \"\\t\" '{{ if (NR!=1){{ print $1 }}}}' {input.genes} > {output}"""

rule genes_ripgrep_geneid:
    input:
        geneid="results/2019-06-26/dge/rg/{cell}_de_geneid.txt",
        refseq="data/2018-11-09/hg19/genes.gtf",
    output:
        "results/2019-06-26/dge/rg/{cell}_de_ripgrep.gtf",
    conda:
        "../envs/ripgrep.yaml"
    threads: 4
    shell:
        """rg --dfa-size-limit 2G -w -f {input.geneid} {input.refseq} > {output}"""

rule genes_ripgrep_exon:
    input:
        rg="results/2019-06-26/dge/rg/{cell}_de_ripgrep.gtf",
    output:
        "results/2019-06-26/dge/rg/{cell}_de_genes.gtf",
    conda:
        "../envs/gawk.yaml"
    threads: 2
    shell:
        """awk '$3 ~ /^exon$/ ' {input} > {output}"""

rule genes_NOIseq:
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

rule groHMM:
    input:
        "results/2019-06-03/eRNA/counts/{sample}_merged.txt"
    output:
        "results/2019-06-26/dge/grohmm/{sample}_grohmm.txt"
    conda:
        "../../envs/grohmm.yaml"
    threads: 4
    script:
        "../../scripts/grohmm.R"
