rule IMR_hg19_bowtie2:
    input:
        sample=["data/2018-06-23/IMR/{unit}.fastq"]
    output:
        "results/2018-12-01/IMR/{unit}.bam"
    log:
        "logs/IMR/bowtie2/{unit}.log"
    params:
        index="data/2018-11-27/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"
