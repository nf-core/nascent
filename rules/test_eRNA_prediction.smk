rule test_eRNA_vs_Peng:
    input:
        edmund="results/2018-12-02/{genome}/{cell}_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
    output:
        report("results/2018-11-10/test/{genome}/{cell}_eRNA_overlaps.bed", category="eRNA Prediction")
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
        report("results/2018-11-10/test/hg19_eRNA_overlaps.bed", category="eRNA Prediction")
    log:
        "logs/hg18/test_GM18_liftOver_eRNAs.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_hg19_vs_hg18:
    input:
        hg19="results/2018-12-02/hg19/{cell}_eRNA.bed",
        hg18="results/2018-12-02/hg18/{cell}_eRNA.bed",
    output:
        report("results/2018-11-10/test/{cell}_hg19_vs_hg18_eRNA.bed", category="eRNA Prediction")
    log:
        "logs/hg19/test_{cell}_eRNAs_hg19_vs_hg18.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.hg19} -b {input.hg18} \
        -sorted -u > {output} 2> {log}"

rule test_eRNA_vs_liftOver:
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        eRNA="results/2018-12-02/{genome}/{cell}_eRNA.bed",
    output:
        report("results/2018-11-10/test/{genome}/{cell}_eRNA_vs_liftOver.bed", category="eRNA Prediction")
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
        report("results/2018-11-10/test/IMR_eRNA_vs_GM19.bed", category="eRNA Prediction")
    log:
        "logs/hg19/test_IMR_vs_GM.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.GM19} \
        -sorted -u > {output} 2> {log}"
