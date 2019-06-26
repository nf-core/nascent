rule GM18_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/GM18/{unit}.bam",unit=GM_SAMPLES)
    output:
        directory("results/2018-11-07/GM18_meta_tagDir")
    singularity:
        "docker://emiller88/homer:latest"
    shell:
        "makeTagDirectory {output} -genome hg18 -checkGC {input}"

rule GM18_meta_findPeaks:
    input:
        "results/2018-11-07/GM18_meta_tagDir/"
    output:
        "results/2018-11-07/GM18_meta_groseq_peak.gtf"
    singularity:
        "docker://emiller88/homer:latest"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule GM18_meta_pos2bed:
    input:
        "results/2018-11-07/GM18_meta_groseq_peak.gtf",
    output:
        "results/2018-11-07/GM18_meta_groseq_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"
