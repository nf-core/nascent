rule IMR_fastqc:
    input:
        "data/2018-06-23/IMR/{sample}.fastq"
    output:
        html="results/2018-10-03/fastqc/IMR/{sample}.html",
        zip="results/2018-10-03/fastqc/IMR/{sample}.zip"
    params: ""
    log:
        "logs/IMR/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
