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

rule GM19_eRNAdiffExpression:
    input:
        "results/2019-02-05/GM/GM19_eRNA_outputannotation.txt"
    output:
        "results/2019-02-05/GM/diffOutput.genes.txt"
    shell:
        "getDiffExpression.pl {input} {GM_SAMPLES} -edgeR -simpleNorm -dispersion 0.05 > {output}"

# Lacks Normalization
rule GM19_diff_Peaks:
    input:
        eRNAPeak="results/2019-01-31/GM19_eRNA_peak.gtf",
        sampleTags="results/2019-01-28/GM/{unit}_tagDir/",
        GM0hTags="results/2019-01-28/GM/GM0h_tagDir/",
    output:
        "results/2019-02-01/GM_annotations/{unit}_diffOutput.tsv"
    params:
    shell:
        "getDifferentialPeaks {input.eRNAPeak} {input.sampleTags} {input.GM0hTags} > {output}"
