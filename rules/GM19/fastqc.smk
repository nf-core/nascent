rule GM18_fastqc:
    input:
        "data/2018-06-23/GM/{sample}.fastq"
    output:
        html="results/2018-10-03/fastqc/GM/{sample}.html",
        zip="results/2018-10-03/fastqc/GM/{sample}.zip"
    params: ""
    log:
        "logs/GM19/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
