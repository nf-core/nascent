rule IMR_RemoveGenes:
    input:
        IMR="results/2018-11-07/IMR_meta_groseq_peak.bed",
        refseq="data/2018-11-09/hg19/genes.gtf",
    output:
        "results/2018-11-09/IMR_groseq_noGenes.bed"
    log:
        "logs/IMR/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"
