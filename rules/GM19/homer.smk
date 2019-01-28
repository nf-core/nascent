rule GM19_meta_makeTagDirectory:
    input:
        expand("results/2018-11-27/GM19/{unit}.bam",unit=GM_SAMPLES)
    output:
        "results/2018-11-28/GM19_meta_tagDir/"
    # conda:
    #     "../../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM19_meta_findPeaks:
    input:
        "results/2018-11-28/GM19_meta_tagDir/"
    output:
        "results/2018-11-28/GM19_meta_groseq_peak.gtf"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule GM19_meta_pos2bed:
    input:
        "results/2018-11-28/GM19_meta_groseq_peak.gtf"
    output:
        "results/2018-11-28/GM19_meta_groseq_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"
