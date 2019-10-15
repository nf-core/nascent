# FIXME uses unmerged regions of eRNAs
rule eRNA_saf_cell:
    input:
        "results/2018-12-02/{genome}/{cell}_eRNA.bed"
    output:
        "results/2019-06-07/{genome}/{cell}_eRNA.saf",
    conda:
        "../../envs/pandas.yaml"
    log:
        "logs/{genome}/{cell}/eRNA_saf.log"
    script:
        "../../scripts/bed2saf.py"

rule eRNA_saf_viral:
    input:
        "results/2019-08-26/eRNA_viral.bed",
    output:
        "results/2019-09-27/eRNA_viral.saf",
    conda:
        "../../envs/pandas.yaml"
    log:
        "logs/eRNA_saf_viral.log"
    script:
        "../../scripts/bed2saf.py"
