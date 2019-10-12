rule test_GM18_vs_Peng:
    input:
        edmund="results/2018-12-02/hg18/GM_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-10/test/hg18_eRNA_overlaps.bed"
    log:
        "logs/GM18/test_GM18_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_GM18_liftOver_vs_Peng:
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-10/test/hg19_eRNA_overlaps.bed"
    log:
        "logs/GM18/test_GM18_liftOver_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"
