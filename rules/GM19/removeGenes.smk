# Add 1 KB to the front of the genes and 10 KB to the end
# | 1 KB | GENE |      10 KB      |

rule hg19_slopRefSeq:
    input:
        refSeq="data/2018-11-09/hg19/genes.bed",
        chromLen="data/2018-06-24/hg19/chrom.sizes"
    output:
        "data/2018-11-09/hg19/hg19_slop_refseq.bed"
    conda:
        "../../envs/bedtools.yaml"
    shell:
         "slopBed -i {input.refSeq} \
         -g {input.chromLen} \
         -l 1000 -r 10000 > {output}"

rule hg19_fixBEDcoordinates:
    input:
        "data/2018-11-09/hg19/hg19_slop_refseq.bed"
    output:
        "data/2018-11-09/hg19/hg19_slop_refseq.sorted.bed"
    log:
        "logs/GM19/RemoveGenes.log"
    conda:
        "../../envs/bedops.yaml"
    shell:
         "awk '{{ if ($2 > $3) {{ t = $2; $2 = $3; $3 = t; }} \
         else if ($2 == $3) {{ $3 += 1; }} print $0; }}' OFS='\\t' \
         {input} | sort-bed - > {output}"

# Take the slopped genes and remove them from the Peaks file
rule GM19_RemoveGenes:
    input:
        GM="results/2018-11-07/GM19_meta_groseq_peak.bed",
        refseq="data/2018-11-09/hg19/hg19_slop_refseq.sorted.bed",
    output:
        "results/2018-11-09/GM19_meta_groseq_noGenes.bed"
    log:
        "logs/GM19/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.GM} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"
