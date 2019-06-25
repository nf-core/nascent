rule GM19_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/GM19/{unit}.bam",unit=GM_SAMPLES)
    output:
        directory("results/2018-11-07/GM19_meta_tagDir/")
    singularity:
        "docker://emiller88/homer:latest"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM19_meta_findPeaks:
    input:
        "results/2018-11-07/GM19_meta_tagDir/"
    output:
        "results/2018-11-07/GM19_meta_groseq_peak.gtf"
    singularity:
        "docker://emiller88/homer:latest"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule GM19_meta_pos2bed:
    input:
        "results/2018-11-07/GM19_meta_groseq_peak.gtf"
    output:
        "results/2018-11-07/GM19_meta_groseq_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"

rule GM19_sample_makeTagDirectory:
    input:
        sample=["results/2018-10-04/GM19/{unit}.bam"],
    output:
        directory("results/2019-01-28/GM/{unit}_tagDir/"),
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM19_sample_findPeaks:
    input:
        "results/2019-01-28/GM/{unit}_tagDir/",
    output:
        "results/2019-01-28/GM/{unit}_groseq_peak.gtf"
    singularity:
        "docker://emiller88/homer:latest"
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule GM19_sample_pos2bed:
    input:
        "results/2019-01-28/GM/{unit}_groseq_peak.gtf"
    output:
        "results/2019-01-28/GM/{unit}_groseq_peak.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"
