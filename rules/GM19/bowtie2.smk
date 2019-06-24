rule GM19_bowtie2:
    input:
        sample=["data/2018-06-23/GM/{unit}.fastq"]
    output:
        "results/2018-10-04/GM19/{unit}.bam"
    log:
        "logs/GM19/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/hg19/genome",
        extra=""
    threads: 4
    wrapper:
        "0.35.1/bio/bowtie2/align"
