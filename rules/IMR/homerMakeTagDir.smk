rule IMR_hg19_All_makeTagDirectory:
    input:
        expand("results/2018-12-01/IMR/{unit}.bam",unit=IMR_SAMPLES)
    output:
        "results/2018-12-02/All_together/"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"
