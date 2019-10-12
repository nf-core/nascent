rule test_eRNA_vs_Peng:
    input:
        edmund="results/2018-12-02/{genome}/{cell}_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
    output:
        "results/2018-11-10/test/{genome}/{cell}_eRNA_overlaps.bed"
    log:
        "logs/{genome}/test_{cell}_eRNAs.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_GM18_liftOver_vs_Peng:
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
    output:
        "results/2018-11-10/test/hg19_eRNA_overlaps.bed"
    log:
        "logs/hg18/test_GM18_liftOver_eRNAs.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_GM19_vs_GM18:
    input:
        GM19="results/2018-12-02/hg19/GM_eRNA.bed",
        GM18="results/2018-12-02/hg18/GM_eRNA.bed",
    output:
        "results/2018-11-10/test/GM19_eRNA_vs_GM18.bed"
    log:
        "logs/hg19/test_GM19_eRNAs_GM18.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.GM19} -b {input.GM18} \
        -sorted -u > {output} 2> {log}"

rule test_eRNA_vs_liftOver:
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        eRNA="results/2018-12-02/{genome}/{cell}_eRNA.bed",
    output:
        "results/2018-11-10/test/{genome}/{cell}_eRNA_vs_liftOver.bed"
    log:
        "logs/{genome}/test_{cell}_eRNAs_liftover.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.eRNA} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_IMR_vs_GM19:
    input:
        IMR="results/2018-12-02/hg19/IMR_eRNA.bed",
        GM19="results/2018-12-02/hg19/GM_eRNA.bed",
    output:
        "results/2018-11-10/test/IMR_eRNA_vs_GM19.bed"
    log:
        "logs/hg19/test_IMR_vs_GM.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.GM19} \
        -sorted -u > {output} 2> {log}"
