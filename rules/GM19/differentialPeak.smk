# TODO Plug this in
rule GM19_merge_sample_peaks:
    input:
        expand("results/2019-01-28/GM/{unit}_groseq_peak.gtf",unit=GM_SAMPLES),
    output:
        "results/2019-01-31/GM/GM19_merged_peaks.gtf"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    params:
        maxDistance="100"
    shell:
        "mergePeaks -d {params.maxDistance} {input} > {output}"

rule GM19_eRNA_makeTagDirectory:
    input:
        "results/2018-11-29/GM19_eRNA.bed",
    output:
        "results/2019-01-31/GM19_eRNA_tagDir/",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM19_sample_DiffPeaksReplicates:
    input:
        sampleTags=expand("results/2019-01-28/GM/{unit}_tagDir/",unit=GM_SAMPLES),
        eRNATags="results/2019-01-31/GM19_eRNA_tagDir/",
    output:
        "results/2019-01-31/GM19_eRNA_diffPeaks.txt",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer_R.yaml"
    params:
        genome="hg19",
    shell:
        "getDifferentialPeaksReplicates.pl -t {input.eRNATags} \ "
        "-i {input.sampleTags} -genome {params.genome} > {output}"
