rule IMR_hg19_bowtie2:
    input:
        sample=["data/2018-06-23/IMR/{unit}.fastq"]
    output:
        "results/2018-10-04/IMR/{unit}.bam"
    log:
        "logs/IMR/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/hg19/genome",
        extra=""
    threads: 4
    wrapper:
        "0.35.1/bio/bowtie2/align"

rule IMR_bam_bed:
    input:
        "results/2018-10-04/IMR/{unit}.bam"
    output:
        "results/2019-10-06/IMR/{unit}.bed"
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/IMR/bamtobed/{unit}.log"
    threads: 2
    shell:
        "bamToBed -i {input} | sort -k1,1 -k2,2n > {output}"

rule bed_check_L2:
    """
    Check for the presence of L2
    chr9 21096515 21098279
    """
    input:
        l2="data/2019-10-01/L2.bed",
        sample=["results/2019-10-06/IMR/{unit}.bed"],
    output:
        "results/2019-10-06/l2/{unit}_bam_L2.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/genes/{unit}_bam_check_L2.log"
    threads: 2
    shell:
        "bedtools closest -d -a {input.l2} -b {input.sample} > {output}"
