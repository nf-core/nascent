# Add 1 KB to the front of the genes and 10 KB to the end
# | 1 KB | GENE |      10 KB      |

rule hg19_slopRefSeq:
    input:
        refSeq="data/2018-11-28/genes.bed",
        chromLen="data/2018-11-27/genome.fa.fai"
    output:
        "results/2018-11-29/sloprefseqhg19.bed"
    conda:
        "../../envs/bedtools.yaml"
    shell:
         "slopBed -i {input.refSeq} \
         -g {input.chromLen} \
         -l 1000 -r 10000 > {output}"

rule hg19_fixBEDcoordinates:
    input:
        "results/2018-11-29/sloprefseqhg19.bed"
    output:
        "results/2018-11-29/sloprefseqhg19.sorted.bed"
    log:
        "logs/RemoveGenes.log"
    conda:
        "../../envs/bedops.yaml"
    shell:
         "awk '{{ if ($2 > $3) {{ t = $2; $2 = $3; $3 = t; }} \
         else if ($2 == $3) {{ $3 += 1; }} print $0; }}' OFS='\\t' \
         {input} | sort-bed - > {output}"

# Take the slopped genes and remove them from the Peaks file
rule hg19_RemoveGenes:
    input:
        GM="results/2018-11-28/GM19_meta_groseq_peak.bed",
        refseq="results/2018-11-29/sloprefseqhg19.sorted.bed",
    output:
        "results/2018-11-29/GM19_meta_groseq_noGenes.bed"
    log:
        "logs/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.GM} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"
