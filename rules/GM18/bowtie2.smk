def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("data/2018-06-23/{unit}.fastq",
                       **wildcards)
    # single end sample
    return "data/2018-06-23/{unit}.fastq".format(**wildcards)


rule bowtie2:
    input:
        sample=["data/2018-06-23/GM/{unit}.fastq"]
    output:
        "results/2018-10-04/GM/{unit}.bam"
    log:
        "logs/bowtie2/{unit}.log"
    params:
        index="data/2018-06-24/genome",
        extra=""
    threads: 4
    wrapper:
        "0.27.1/bio/bowtie2/align"
