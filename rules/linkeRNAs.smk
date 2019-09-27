rule eRNA_link_genes:
    input:
        # FIXME Just refseq right now
        dges="results/2019-06-26/dge/rg/{cell}_de_genes.gtf",
        eRNA="results/2018-12-02/{cell}_eRNA.bed",
    output:
        "results/2019-08-26/{cell}_link_eRNA.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
        window="200000"
    threads: 2
    shell:
        "bedtools window -u -w {params.window} -a {input.eRNA} -b {input.dges} > {output}"

rule eRNA_link_overlap:
    input:
        GM="results/2019-08-26/GM19_link_eRNA.bed",
        IMR="results/2019-08-26/IMR_link_eRNA.bed",
    output:
        "results/2019-08-26/overlap_link_eRNA.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
    threads: 3
    shell:
        "bedtools intersect -wo -a {input.GM} -b {input.IMR} > {output}"
