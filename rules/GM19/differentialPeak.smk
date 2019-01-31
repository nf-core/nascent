rule GM19_eRNA_makeTagDirectory:
    input:
        "results/2018-11-29/GM19_eRNA.bed",
    output:
        "results/2018-01-30/GM19_eRNA_tagDir/",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule GM19_sample_DiffPeaksReplicates:
    input:
        sampleTags=expand("results/2019-01-28/GM/{unit}_tagDir/",unit=GM_SAMPLES),
        eRNATags="results/2018-01-30/GM19_eRNA_tagDir/",
    output:
        "results/2018-01-30/GM19_eRNA_diffPeaks.txt",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer_R.yaml"
    params:
        genome="hg19",
    shell:
        "getDifferentialPeaksReplicates.pl -t {input.eRNATags} \ "
        "-i {input.sampleTags} -genome hg19 > {output}"
