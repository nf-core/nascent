rule IMR_hg19_meta_makeTagDirectory:
    input:
        expand("results/2018-12-01/IMR/{unit}.bam",unit=IMR_SAMPLES)
    output:
        directory("results/2018-12-02/IMR_meta_tagDir/")
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule IMR_hg19_meta_findPeaks:
    input:
        "results/2018-12-02/IMR_meta_tagDir/",
    output:
        "results/2018-12-02/IMR_meta_groseq_peak.gtf",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule IMR_hg19_meta_pos2bed:
    input:
        "results/2018-12-02/IMR_meta_groseq_peak.gtf",
    output:
        "results/2018-12-02/IMR_meta_groseq_peak.bed",
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"

rule IMR_hg19_sample_makeTagDirectory:
    input:
        sample=["results/2018-12-01/IMR/{unit}.bam"],
    output:
        "results/2019-01-28/IMR/{unit}_tagDir/",
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input.sample}"

rule IMR_hg19_sample_findPeaks:
    input:
        "results/2019-01-28/IMR/{unit}_tagDir/",
    output:
        "results/2019-01-28/IMR/{unit}_groseq_peak.gtf"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule IMR_hg19_sample_pos2bed:
    input:
        "results/2019-01-28/IMR/{unit}_groseq_peak.gtf"
    output:
        "results/2019-01-28/IMR/{unit}_groseq_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"
