rule test_hg18_eRNAs:
    input:
        edmund="results/2018-11-10/eRNA_GM_hg18.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-10/hg18_eRNA_overlaps.bed"
    log:
        "logs/test_hg18_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_hg19_eRNAs:
    input:
        edmund="results/2018-11-25/eRNA_GM_liftover_hg19.sorted.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.sorted.bed",
    output:
        "results/2018-11-25/hg19_eRNA_overlaps.bed"
    log:
        "logs/test_hg19_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"
