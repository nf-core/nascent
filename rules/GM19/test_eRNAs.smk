rule test_GM19_vs_Peng:
    input:
        edmund="results/2018-12-02/hg19/GM_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-10/test/GM19_eRNA_vs_Peng.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_GM19_vs_GM18:
    input:
        GM19="results/2018-12-02/hg19/GM_eRNA.bed",
        GM18="results/2018-12-02/hg18/GM_eRNA.bed",
    output:
        "results/2018-11-10/test/GM19_eRNA_vs_GM18.bed"
    log:
        "logs/GM19/test_GM19_eRNAs_GM18.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.GM19} -b {input.GM18} \
        -sorted -u > {output} 2> {log}"

rule test_GM19_vs_liftOver:
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        GM19="results/2018-12-02/hg19/GM_eRNA.bed",
    output:
        "results/2018-11-10/test/GM19_eRNA_vs_liftOver.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.GM19} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"
