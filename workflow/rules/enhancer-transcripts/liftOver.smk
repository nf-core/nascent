rule GM18_liftOver_eRNAs:
    """
    Calls liftover on GM eRNAs that were predictied using hg18 to hg19
    """
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
    """
    Sorts the lifted over eRNAs
    """
    input:
        "results/2018-11-10/eRNA_GM_liftover_hg19.bed",
    output:
        "results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "sortBed -i {input} > {output}"
