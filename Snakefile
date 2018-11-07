# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


# configfile: "config.yaml"


rule all:
    input:
        html=expand("results/2018-10-03/fastqc/{sample}.html",samples=SAMPLES),
        zip=expand("results/2018-10-03/fastqc/{sample}.zip",samples=SAMPLES),

include: "rules/fastqc.smk"
