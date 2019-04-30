rule GM19_bowtie2:
    input:
        sample=["data/2018-06-23/GM/{unit}.fastq"]
    output:
        "results/2018-11-27/GM19/{unit}.bam"
    log:
        "logs/GM19/bowtie2/{unit}.log"
    params:
        index="data/2018-11-27/genome",
        extra=""
    threads: 8
    conda: "../../envs/bowtie2.yaml"
    wrapper:
        "0.27.1/bio/bowtie2/align"
