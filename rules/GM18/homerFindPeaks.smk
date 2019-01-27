rule findPeaks:
    input:
        "results/2018-11-07/All_together/"
    output:
        "results/2018-11-07/groseq_GM_peak.gtf"
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"
