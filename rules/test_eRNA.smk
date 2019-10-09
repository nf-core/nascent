rule test_eRNA_check_L2:
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

rule test_homer_check_L2:
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

rule test_report_L2:
    input:
        expand("results/2019-08-26/{cell}_L2.bed", cell=["GM19","IMR"]),
        expand("results/2019-10-01/{cell}_genes_L2.bed", cell=["GM19","IMR"]),
    output:
        report("results/2019-10-01/L2.csv", caption="../report/L2.rst", category="eRNA")
    log:
        "logs/report/L2.log"
    threads: 4
    shell:
        "cat {input} > {output}"
