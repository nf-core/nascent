rule GM18_Meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/GM18/{unit}.bam",unit=GM_SAMPLES)
    output:
        "results/2018-11-07/GM18_meta_tagDir/"
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg18 -checkGC {input}"
