rule bowtie2:
    input:
        sample=["data/2018-06-23/{unit}.fastq"]
    output:
        "results/2018-10-04/{unit}.bam"
    log:
        "logs/bowtie2/{unit}.log"
    params:
        index=config["ref"]["index"],
        extra=""
    threads: 4
    wrapper:
        "0.38.0/bio/bowtie2/align"
