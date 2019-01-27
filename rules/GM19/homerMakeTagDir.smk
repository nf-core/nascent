rule GM_hg19_All_makeTagDirectory:
    input:
        expand("results/2018-11-27/GM/{unit}.bam",unit=GM_SAMPLES)
    output:
        "results/2018-11-28/All_together/"
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"
