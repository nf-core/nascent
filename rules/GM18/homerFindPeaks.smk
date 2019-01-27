rule GM18_findPeaks:
    input:
        "results/2018-11-07/GM18_meta_tagDir/"
    output:
        "results/2018-11-07/GM18_meta_groseq_peak.gtf"
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"
