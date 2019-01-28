rule GM_hg19_All_makeTagDirectory:
    input:
        expand("results/2018-11-27/GM19/{unit}.bam",unit=GM_SAMPLES)
    output:
        "results/2018-11-28/GM19_meta_tagDir/"
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"
