rule IMR_eRNA_peak:
    input:
        "results/2018-12-02/eRNA_IMR_hg19.bed",
    output:
        "results/2019-02-05/IMR_eRNA_peak.gtf",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "bed2pos.pl {input} > {output}"

rule IMR_eRNA_analyzeRepeats:
    input:
        eRNAs="results/2019-02-05/IMR_eRNA_peak.gtf",
        IMR_tags=expand("results/2019-01-28/IMR/{unit}_tagDir/", unit=IMR_SAMPLES),
    output:
        "results/2019-02-05/IMR/IMR_eRNA_outputannotation.txt"
    params:
        genome="hg19",
    shell:
        "analyzeRepeats.pl {input.eRNAs} {params.genome} -count genes -raw -d {input.IMR_tags} > {output}"

rule IMR_eRNAdiffExpression:
    input:
        "results/2019-02-05/IMR/IMR_eRNA_outputannotation.txt",
    output:
        "results/2019-02-05/IMR/IMR_eRNA_diffOutput.txt"
    shell:
        "getDiffExpression.pl {input} {IMR_SAMPLES} -edgeR -simpleNorm -dispersion 0.05 > {output}"

# OVERLAP
rule overlap_eRNA_peak:
    input:
        "results/2018-12-03/IMR_hg19_vs_GM_hg19.bed"
    output:
        "results/2019-02-05/IMR_hg19_vs_GM19.gtf",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "bed2pos.pl {input} > {output}"

rule IMR_overlap_eRNA_analyzeRepeats:
    input:
        eRNAs="results/2019-02-05/IMR_hg19_vs_GM19.gtf",
        IMR_tags=expand("results/2019-01-28/IMR/{unit}_tagDir/", unit=IMR_SAMPLES),
    output:
        "results/2019-02-05/IMR/IMR_overlap_eRNA_outputannotation.txt"
    params:
        genome="hg19",
    shell:
        "analyzeRepeats.pl {input.eRNAs} {params.genome} -count genes -raw -d {input.IMR_tags} > {output}"

rule IMR_overlap_eRNAdiffExpression:
    input:
        "results/2019-02-05/IMR/IMR_overlap_eRNA_outputannotation.txt"
    output:
        "results/2019-02-05/IMR/IMR_overlap_eRNA_diffOutput.txt"
    shell:
        "getDiffExpression.pl {input} {IMR_SAMPLES} -edgeR -simpleNorm -dispersion 0.05 > {output}"
