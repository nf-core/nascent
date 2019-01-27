rule GM18_fastqc:
    input:
        "data/2018-06-23/GM18/{sample}.fastq"
    output:
        html="results/2018-10-03/GM18_fastqc/{sample}.html",
        zip="results/2018-10-03/GM18_fastqc/{sample}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
