rule test_eRNA_check_L2:
    """
    Check for the presence of L2 in eRNAs
    chr9:21096515-21098279
    """
    input:
        eRNA="results/2019-08-26/hg19/{cell}_link_eRNA.bed",
        l2="reference/2019-10-01/L2.bed",
    output:
        "results/2019-08-26/{cell}_L2.bed",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/eRNA/{cell}_check_L2.log"
    threads: 2
    shell:
        "bedtools closest -a {input.l2} -b {input.eRNA} > {output}"

rule test_homer_check_L2:
    """
    Check for the presence of L2 in GRO-Seq transcripts from homer
    chr9:21096515-21098279
    """
    input:
        homer="results/2018-11-07/{genome}/{cell}_meta_transcripts.bed",
        l2="reference/2019-10-01/L2.bed",
    output:
        "results/2019-10-01/{genome}/{cell}_genes_L2.bed",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/{genome}/genes/{cell}_check_L2.log"
    shell:
        "bedtools closest -a {input.l2} -b {input.homer} > {output}"

rule test_report_L2:
    """
    Reports the nearest transcripts to Legacy L2
    FIXME Possibly add L2 over Histones
    """
    input:
        eRNA=expand("results/2019-08-26/{cell}_L2.bed", cell=["GM","IMR"]),
        genes=expand("results/2019-10-01/hg19/{cell}_genes_L2.bed", cell=["GM","IMR"]),
    output:
        report("results/2019-10-01/L2.tsv", caption="../../report/L2.rst", category="Inducible Pairs")
    log:
        "logs/hg19/test/L2.log"
    shell:
        "cat {input.eRNA} {input.genes} > {output}"

rule fig_linked_eRNA_cross_cell:
    """
    Creates a venn diagram of viral linked eRNAs across cell lines
    """
    input:
        IMR="results/2019-08-26/hg19/IMR_link_eRNA.bed",
        GM="results/2019-08-26/hg19/GM_link_eRNA.bed",
        overlap="results/2019-08-26/eRNA_overlap_viral.bed",
    output:
        report("results/2018-10-01/{genome}/eRNA_cross_cell_viral.svg", category="Figures")
    log:
        "logs/{genome}/figure/eRNA_cross_cell.log"
    conda:
        "../../envs/venn.yaml"
    params:
        title="IMR vs GM Linked eRNAs",
    script:
        "../../scripts/venn-smk.py"
