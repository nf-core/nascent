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
        eRNAPeak="results/2019-01-31/GM19_eRNA_peak.gtf",
        count="results/2019-01-31/GM_annotations/{unit}_countTable.tsv",
    output:
        "results/2019-01-31/GM_annotations/{unit}_diffExpression.txt"
    params:
    shell:
        "getDiffExpression.pl {input} eRNA -peaks > {output}"

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
