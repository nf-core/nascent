rule test_hg19_eRNAs_Peng:
    input:
        edmund="results/2018-11-29/GM19_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-30/hg19_eRNA_overlaps_Peng.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_hg19_eRNAs_liftOver:
    input:
        liftOver="results/2018-11-25/eRNA_GM_liftover_hg19.sorted.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-30/hg19_eRNA_overlaps_liftOver.bed"
    log:
        "logs/GM19/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"
