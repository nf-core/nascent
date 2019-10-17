# TODO Generalize these
rule IMR18_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/hg18/{unit}.bam",unit=IMR_SAMPLES)
    output:
        directory("results/2018-11-07/hg18/IMR_meta_tagDir"),
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "makeTagDirectory {output} -genome hg18 -checkGC {input}"

rule GM18_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/hg18/{unit}.bam",unit=GM_SAMPLES)
    output:
        directory("results/2018-11-07/hg18/GM_meta_tagDir")
    singularity:
        "docker://emiller88/homer:latest"
    threads: 4
    shell:
        "makeTagDirectory {output} -genome hg18 -checkGC {input}"

rule IMR_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/hg19/{unit}.bam",unit=IMR_SAMPLES)
    output:
        directory("results/2018-11-07/hg19/IMR_meta_tagDir"),
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM_meta_makeTagDirectory:
    input:
        expand("results/2018-10-04/hg19/{unit}.bam",unit=GM_SAMPLES)
    output:
        directory("results/2018-11-07/hg19/GM_meta_tagDir")
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule meta_makeUCSCfile:
    input:
        "results/2018-11-07/{genome}/{cell}_meta_tagDir"
    output:
        report("results/2018-11-07/{genome}/{cell}_{genome}_ucsc.zip", category="homer")
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "makeUCSCfile {input} -o {output} -strand separate"

rule hg18_meta_findPeaks:
    input:
        "results/2018-11-07/hg18/{cell}_meta_tagDir"
    output:
        "results/2018-11-07/hg18/{cell}_meta_transcripts.txt"
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "findPeaks {input} -o {output} -style groseq"

rule hg19_meta_findPeaks:
    input:
        tagdir="results/2018-11-07/hg19/{cell}_meta_tagDir",
        uniqmap="data/2019-07-26/hg19-50nt-uniqmap",
    output:
        "results/2018-11-07/hg19/{cell}_meta_transcripts.txt"
    singularity:
        "docker://emiller88/homer:latest"
    threads: 2
    shell:
        "findPeaks {input.tagdir} -style groseq -o {output} -uniqmap {input.uniqmap}"

rule homer_meta_pos2bed:
    input:
        "results/2018-11-07/{genome}/{cell}_meta_transcripts.txt"
    output:
        "results/2018-11-07/{genome}/{cell}_meta_transcripts.bed"
    conda:
        "../../envs/homer.yaml"
    threads: 2
    shell:
        "pos2bed.pl {input} | sort -k1,1 -k2,2n - > {output}"

# rule homer_sample_makeTagDirectory:
#     input:
#         sample=["results/2018-10-04/{unit}.bam"],
#     output:
#         directory("results/2019-01-28/{unit}_tagDir"),
#     threads: 2
#     shell:
#         "makeTagDirectory {output} -genome hg19 -checkGC {input}"

# rule homer_sample_findPeaks:
#     input:
#         tagdir="results/2019-01-28/{unit}_tagDir",
#         uniqmap="data/2019-07-26/hg19-50nt-uniqmap",
#     output:
#         "results/2019-01-28/hg19/{unit}_groseq_peak.gtf"
#     singularity:
#         "docker://emiller88/homer:latest"
#     threads: 2
#     shell:
#         "findPeaks {input.tagdir} -o {output} -style groseq -uniqmap {input.uniqmap}"

# rule homer_sample_pos2bed:
#     input:
#         "results/2019-01-28/hg19/{unit}_groseq_peak.gtf"
#     output:
#         "results/2019-01-28/hg19/{unit}_groseq_peak.bed"
#     conda:
#         "../envs/homer.yaml"
#     shell:
#         "pos2bed.pl {input} > {output}"
