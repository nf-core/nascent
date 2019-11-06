rule slopRefSeq:
    """
    Add 1 KB to the front of the genes and 10 KB to the end
    | 1 KB | GENE |      10 KB      |
    """
    input:
        refSeq="data/2018-11-09/{genome}/genes.bed",
        chromLen="data/2019-09-25/{genome}.chrom.sizes",
    output:
        "data/2018-11-09/{genome}/{genome}_slop_refseq.bed"
    conda:
        "../../envs/bedtools.yaml"
    params:
        upstream=config["slopRefSeq"]["upstream"],
        downstream=config["slopRefSeq"]["downstream"],
    shell:
         "sortBed -i {input.refSeq} |\
         slopBed -i - \
         -g {input.chromLen} \
         -l {params.upstream} -r {params.downstream} > {output}"

rule fixBEDcoordinates:
    """
    Removes chrXYZ_random entries in hg18
    Fixes alignments if start is after the end
    Sorts the entries
    """
    input:
        "data/2018-11-09/{genome}/{genome}_slop_refseq.bed"
    output:
        "data/2018-11-09/{genome}/{genome}_slop_refseq.sorted.bed"
    log:
        "logs/{genome}/RemoveGenes.log"
    conda:
        "../../envs/bedops.yaml"
    shell:
         """
         awk -F '\\t' 'length($1) <= 5 {{ print }}' {input} |
         awk '{{ if ($2 > $3) {{ t = $2; $2 = $3; $3 = t; }} \
         else if ($2 == $3) {{ $3 += 1; }} print $0; }}' OFS='\\t' - \
         | sort-bed - > {output} 2> {log}
         """

rule removeGenes:
    """
    Removes the intergenic regions from GRO-Seq Transcripts
    """
    input:
        transcripts="results/2018-11-07/{genome}/{cell}_meta_transcripts.bed",
        refseq="data/2018-11-09/{genome}/{genome}_slop_refseq.sorted.bed",
    output:
        "results/2018-11-09/{genome}/{cell}_meta_transcripts_noGenes.bed"
    log:
        "logs/{genome}/{cell}/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.transcripts} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"

rule sample_removeGenes:
    """
    Removes the intergenic regions from GRO-Seq Transcripts
    """
    input:
        transcripts="results/2018-11-07/{genome}/{sample}_transcripts.bed",
        refseq="data/2018-11-09/{genome}/{genome}_slop_refseq.sorted.bed",
    output:
        "results/2018-11-09/{genome}/{sample}_transcripts_noGenes.bed"
    log:
        "logs/{genome}/{sample}/RemoveGenes.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.transcripts} -b {input.refseq} -v \
        | sort -k1,1 -k2,2n - > {output} 2> {log}"
