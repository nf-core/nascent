rule fastqc:
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

# rule report_fastqc:
#     input:
#         expand("results/2018-10-03/fastqc/{sample}.html",sample=IMR_SAMPLES),
#     output:
#         report("results/2018-10-03/fastqc/{sample}.html", category="Quality Control"),
#     params: ""
#     log:
#         "logs/IMR/fastqc/{sample}.log"
#     shell:
#         "cp {input} {output}"
