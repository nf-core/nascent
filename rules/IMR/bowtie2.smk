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
