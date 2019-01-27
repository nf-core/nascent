rule All_makeTagDirectory:
    input:
        expand("results/2018-10-04/GM/{unit}.bam",unit=GM_SAMPLES)
    output:
        "results/2018-11-07/All_together/"
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg18 -checkGC {input}"
