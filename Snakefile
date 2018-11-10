# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

singularity: "docker://continuumio/miniconda3:4.4.10"
# configfile: "config.yaml"

SAMPLES=["GM0h", "GM30min", "GM1h", "GM2h", "GM4h", "GM6h",  "GM9h", "GM12h", "GM18h", "GM24h", "GM48h", "GM72h",]

rule all:
    input:
        html=expand("results/2018-10-03/fastqc/{sample}.html",sample=SAMPLES),
        zip=expand("results/2018-10-03/fastqc/{sample}.zip",sample=SAMPLES),
        bam=expand("results/2018-10-04/{sample}.bam", sample=SAMPLES),
        bed_peaks="results/2018-11-08/groseq_GM_peak.gtf",

include: "rules/fastqc.smk"
include: "rules/bowtie2.smk"
include: "rules/homerMakeTagDir.smk"
include: "rules/homerFindPeaks.smk"
include: "rules/pos2bed.smk"
