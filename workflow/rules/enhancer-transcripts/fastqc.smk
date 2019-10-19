rule fastqc:
    """
    Creates quality reports of the fastq files and makes a report category
    """
    input:
        "data/2018-06-23/{sample}.fastq"
    output:
        html="results/2018-10-03/fastqc/{sample}.html",
        zip=report("results/2018-10-03/fastqc/{sample}.zip", category="Quality Control"),
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"
