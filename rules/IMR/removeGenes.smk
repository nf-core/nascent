rule IMR_hg19_RemoveGenes:
    input:
        IMR="results/2018-12-02/groseq_IMR_hg19_peak.bed",
        refseq="results/2018-11-29/sloprefseqhg19.sorted.bed",
    output:
        "results/2018-12-02/groseq_IMR_hg19_noGenes.bed"
    log:
        "logs/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"
