rule bowtie2:
    """
    Aligns fastq samples with bowtie2 using reference genome
    TODO use an input function to get the location of fastq files
    """
    input:
        sample=["data/2018-06-23/{unit}.fastq"],
        ref="data/2018-06-24/{genome}/genome.fa",
    output:
        "results/2018-10-04/{genome}/{unit}.bam"
    log:
        "logs/{genome}/bowtie2/{unit}.log"
    params:
        # index=config["ref"]["index"],
        index= lambda w: config["ref"]["{}".format(w.genome)]["index"],
        extra=""
    threads: 4
    wrapper:
        "0.38.0/bio/bowtie2/align"

rule star_index:
    input:
        fasta = "data/2018-10-16/fasta/{genome}/genome.fa",
    output:
        directory("data/2018-10-16/{genome}/star")
    message:
        "Testing STAR index"
    threads:
        16
    params:
        extra = ""
    log:
        "logs/star_index_{genome}.log"
    wrapper:
        "0.39.0/bio/star/index"

rule star_se:
    input:
        fq1="data/2018-06-23/{sample}.fastq",
        starthing="data/2018-10-16/{genome}/star",
    output:
        "results/2018-10-04/{genome}/{sample}/Aligned.out.bam"
    log:
        "logs/{genome}/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index= lambda w: config["ref"]["{}".format(w.genome)]["starIndex"],
        # optional parameters
        extra=""
    threads: 16
    wrapper:
        "0.39.0/bio/star/align"
