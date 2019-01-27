rule GM_fastqc:
    input:
        "data/2018-06-23/GM/{sample}.fastq"
    output:
        html="results/2018-10-03/fastqc/{sample}.html",
        zip="results/2018-10-03/fastqc/{sample}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
