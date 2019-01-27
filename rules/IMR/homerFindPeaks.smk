rule IMR_hg19_findPeaks:
    input:
        "results/2018-12-02/All_together/"
    output:
        "results/2018-12-02/groseq_IMR_hg19_peak.gtf"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule IMR_hg19_pos2bed:
    input:
        "results/2018-12-02/groseq_IMR_hg19_peak.gtf"
    output:
        "results/2018-12-02/groseq_IMR_hg19_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"
