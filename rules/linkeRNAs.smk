rule eRNA_link_genes:
    input:
        # FIXME Just refseq right now
        dges="results/2019-06-26/dge/DE/{cell}_de_genes.gtf",
        eRNA="results/2018-12-02/{cell}_eRNA.bed",
    output:
        "results/2019-08-26/{cell}_link_eRNA.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
        window="200000"
    shell:
        "bedtools window -u -w {params.window} -a {input.dges} -b {input.eRNA} > {output}"
