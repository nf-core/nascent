rule findPeaks:
    input:
        "results/2018-11-07/All_together/"
    output:
        "results/2018-11-07/GroseqIMR90peak.gtf"
    shell:
        "findPeaks {input} -o {output} -style groseq"
