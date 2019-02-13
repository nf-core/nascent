rule GM19_eRNA_analyzeRepeats:
    input:
        eRNAs="results/2019-01-31/GM19_eRNA_peak.gtf",
        GM_tags=expand("results/2019-01-28/GM/{unit}_tagDir/", unit=GM_SAMPLES),
    output:
        "results/2019-02-05/GM/GM19_eRNA_outputannotation.txt"
    params:
        genome="hg19",
    shell:
        "analyzeRepeats.pl {input.eRNAs} {params.genome} -count genes -raw -d {input.GM_tags} > {output}"

rule GM19_eRNAdiff_Expression:
    input:
        "results/2019-02-05/GM/GM19_eRNA_outputannotation.txt"
    output:
        "results/2019-02-05/GM/diffOutput.genes.txt"
    shell:
        "getDiffExpression.pl {input} {GM_SAMPLES} -edgeR -simpleNorm -dispersion 0.05 > {output}"

# Overlap
rule GM19_overlap_eRNA_peak:
    input:
        "results/2019-02-06/GM_hg19_vs_IMR_hg19.bed",
    output:
        "results/2019-02-06/GM_hg19_vs_IMR_hg19.gtf",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "bed2pos.pl {input} > {output}"
# Comes from IMR/differentialPeak.smk
rule GM19_overlap_eRNA_analyzeRepeats:
    input:
        eRNAs="results/2019-02-06/GM_hg19_vs_IMR_hg19.gtf",
        GM_tags=expand("results/2019-01-28/GM/{unit}_tagDir/", unit=GM_SAMPLES),
    output:
        "results/2019-02-05/GM/GM_overlap_eRNA_outputannotation.txt"
    params:
        genome="hg19",
    shell:
        "analyzeRepeats.pl {input.eRNAs} {params.genome} -count genes -raw -d {input.GM_tags} > {output}"

rule GM19_overlap_eRNAdiffExpression:
    input:
        "results/2019-02-05/GM/GM_overlap_eRNA_outputannotation.txt"
    output:
        "results/2019-02-05/GM/GM_overlap_eRNA_diffOutput.txt"
    shell:
        "getDiffExpression.pl {input} {GM_SAMPLES} -edgeR -simpleNorm -dispersion 0.05 > {output}"
