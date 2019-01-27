rule test_IMR90_hg19_vs_Peng:
    input:
        IMR="results/2018-12-02/eRNA_IMR_hg19.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-12-03/IMR_hg19_vs_Peng.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.peng} \
        -sorted -u > {output} 2> {log}"

rule test_IMR_hg19_vs_liftOver:
    input:
        GM_liftOver="results/2018-11-29/eRNA_GM_hg19.bed",
        IMR_hg19eRNAs="results/2018-12-02/eRNA_IMR_hg19.bed",
    output:
        "results/2018-12-03/IMR_hg19_vs_liftOver.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR_hg19eRNAs} -b {input.GM_liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_IMR_hg19_vs_GM_hg19:
    input:
        IMR_hg19eRNAs="results/2018-12-02/eRNA_IMR_hg19.bed",
        GM_hg19eRNAs="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-12-03/IMR_hg19_vs_GM_hg19.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR_hg19eRNAs} -b {input.GM_hg19eRNAs} \
        -sorted -u > {output} 2> {log}"
