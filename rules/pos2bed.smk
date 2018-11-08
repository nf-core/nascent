rule pos2bed:
    input:
        "results/2018-11-07/groseq_GM_peak.gtf"
    output:
        "results/2018-11-08/groseq_GM_peak.bed"
    run:
        shell('pos2bed.pl {input} > {output}')
