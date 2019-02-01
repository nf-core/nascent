rule GM19_annotatePeaks_sample:
    input:
        peakFile="results/2019-01-28/GM/{unit}_groseq_peak.gtf",
        eRNAPeak="results/2019-01-31/GM19_eRNA_peak.gtf",
    output:
        "results/2019-01-31/GM_annotations/{unit}_countTable.tsv"
    params:
        genome="hg19",
    shell:
        "annotatePeaks.pl {input.peakFile} {input.eRNAPeak} {params.genome} -raw > {output}"

rule GM19_diffExpression:
    input:
        "results/2019-01-31/GM_annotations/{unit}_countTable.tsv",
    output:
        expand("results/2019-01-31/GM_annotations/{unit}_diffOutput.txt",unit=GM_SAMPLES),
    params:
    shell:
        "getDiffExpression.pl {input} {unit} eRNA -peaks > {output}"

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
