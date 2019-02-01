    input:
    output:
    params:
    shell:

    input:
    output:
    shell:

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
