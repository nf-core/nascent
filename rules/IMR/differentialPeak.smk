rule IMR_eRNA_makeTagDirectory:
    input:
        "results/2018-12-02/eRNA_IMR_hg19.bed"
    output:
        "results/2018-01-30/IMR_eRNA_tagDir/",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule IMR_sample_DiffPeaksReplicates:
    input:
        sampleTags=expand("results/2018-12-01/IMR/{unit}.bam",unit=IMR_SAMPLES),
        eRNATags="results/2018-01-30/IMR_eRNA_tagDir/",
    output:
        "results/2018-01-30/IMR_eRNA_diffPeaks.txt",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer_R.yaml"
    params:
        genome="hg19",
    shell:
        "getDifferentialPeaksReplicates.pl -t {input.eRNATags} \ "
        "-i {input.sampleTags} -genome hg19 > {output}"
