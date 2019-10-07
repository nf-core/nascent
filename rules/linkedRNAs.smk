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

rule eRNA_link_merge:
    input:
        GM="results/2019-08-26/GM19_link_eRNA.bed",
        IMR="results/2019-08-26/IMR_link_eRNA.bed",
    output:
        "results/2019-08-26/eRNA_merged.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
        # dis="-d 2000",
        col="-c 4,5,6 -o distinct",
    threads: 3
    shell:
        "cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin -s {params.col} > {output}"

rule eRNA_check_L2:
    """
    Check for the presence of L2
    chr9 21096515 21098279
    """
    input:
        eRNA="results/2019-08-26/{cell}_link_eRNA.bed",
        # eRNA="results/2018-12-02/{cell}_eRNA.bed",
        l2="data/2019-10-01/L2.bed",
    output:
        "results/2019-08-26/{cell}_L2.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/eRNA/{cell}_check_L2.log"
    threads: 2
    shell:
        # TODO Make a report
        "bedtools closest -a {input.l2} -b {input.eRNA} > {output}"

rule homer_check_L2:
    """
    Check for the presence of L2
    chr9 21096515 21098279
    """
    input:
        homer="results/2018-11-07/{cell}_meta_groseq_peak.bed",
        l2="data/2019-10-01/L2.bed",
    output:
        "results/2019-10-01/{cell}_genes_L2.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/genes/{cell}_check_L2.log"
    threads: 2
    shell:
        "bedtools closest -a {input.l2} -b {input.homer} > {output}"

rule report_L2:
    input:
        expand("results/2019-08-26/{cell}_L2.bed", cell=["GM19","IMR"]),
        expand("results/2019-10-01/{cell}_genes_L2.bed", cell=["GM19","IMR"]),
    output:
        report("L2.tsv", caption="report/L2.rst", category="eRNA")
    log:
        "logs/report/L2.log"
    threads: 4
    shell:
        "cat {input} > {output}"

rule eRNA_link_overlap:
    input:
        GM="results/2019-08-26/GM19_link_eRNA.bed",
        IMR="results/2019-08-26/IMR_link_eRNA.bed",
    output:
        "results/2019-08-26/eRNA_overlap_viral.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
    threads: 3
    shell:
        "bedtools intersect -wo -a {input.GM} -b {input.IMR} > {output}"

rule eRNA_link_area:
    input:
        merge="results/2019-08-26/eRNA_merged.bed",
        ov="results/2019-08-26/eRNA_overlap_viral.bed",
    output:
        "results/2019-08-26/eRNA_viral.bed",
    conda:
        "../envs/bedtools.yaml"
    params:
    threads: 3
    shell:
        "bedtools intersect -wa -a {input.merge} -b {input.ov} > {output}"
