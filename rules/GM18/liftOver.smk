rule GM18_liftOver_eRNAs:
    input:
        eRNAs="results/2018-12-02/hg18/GM_eRNA.bed",
        mapChain="data/2018-11-10/hg18ToHg19.over.chain"
    output:
        "results/2018-11-10/eRNA_GM_liftover_hg19.bed",
    log:
        "logs/GM18/liftover.log"
    conda:
        "../../envs/liftover.yaml"
    shell:
        "liftOver {input.eRNAs} {input.mapChain} {output} unmapped 2> {log}"

rule GM18_sort_liftOver_Bed:
    input:
        "results/2018-11-10/eRNA_GM_liftover_hg19.bed",
    output:
        "results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    threads: 2
    shell:
        "sortBed -i {input} > {output}"
