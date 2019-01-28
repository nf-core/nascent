rule IMR_fastqc:
    input:
        "data/2018-06-23/IMR/{sample}.fastq"
    output:
        html="results/2018-12-01/fastqc/{sample}.html",
        zip="results/2018-12-01/fastqc/{sample}.zip"
    params: ""
    log:
        "logs/IMR/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
