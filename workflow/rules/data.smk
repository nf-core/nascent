# Datasets TODO

# TODO move path to config
rule GM_download:
    output:
        expand("data/2018-06-23/{unit}.fastq",unit=GM_SAMPLES),

# TODO move path to config
rule IMR_download:
    output:
        expand("data/2018-06-23/{unit}.fastq",unit=IMR_SAMPLES),


############
# Histones #
############

# TODO move path to config
# https://www.encodeproject.org/experiments/ENCSR000AKF/
GM_h3k4me1_gz="https://www.encodeproject.org/files/ENCFF000ASM/@@download/ENCFF000ASM.fastq.gz"
# https://www.encodeproject.org/experiments/ENCSR000AKC/
GM_h3k27ac_gz="https://www.encodeproject.org/files/ENCFF000ASU/@@download/ENCFF000ASU.fastq.gz"

rule GM_H3K4me1_fastq:
    output:
        "data/2018-11-13/GM_H3K4me1.fastq",
    shell:
        "curl -Ls {GM_h3k4me1_gz} | gunzip > {output}"

rule GM_H3K27ac_fastq:
    output:
        "data/2018-11-13/GM_H3K27ac.fastq",
    shell:
        "curl -Ls {GM_h3k27ac_gz} | gunzip > {output}"

# TODO Figure out the agregation
# TODO move path to config
# https://www.encodeproject.org/experiments/ENCSR831JSP/
IMR_h3k4me1_gz="https://www.encodeproject.org/files/ENCFF123RFO/@@download/ENCFF123RFO.fastq.gz"
# https://www.encodeproject.org/experiments/ENCSR002YRE/
IMR_h3k27ac_gz="https://www.encodeproject.org/files/ENCFF200XEH/@@download/ENCFF200XEH.fastq.gz"

rule IMR_H3K4me1_fastq:
    output:
        "data/2018-11-13/IMR_H3K4me1.fastq",
    shell:
        "curl -Ls {IMR_h3k4me1_gz} | gunzip > {output}"

rule IMR_H3K27ac_fastq:
    output:
        "data/2018-11-13/IMR_H3K27ac.fastq",
    shell:
        "curl -Ls {IMR_h3k27ac_gz} | gunzip > {output}"

# liftOver
# TODO move path to config
hg18mapChain="http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"

rule download_mapChain:
    output:
        "data/2018-11-10/hg18ToHg19.over.chain"
    shell:
        "curl -Ls {hg18mapChain} | gunzip > {output}"

#########
# Homer #
#########

rule download_hg19_uniqmap:
    output:
        "data/2019-07-26/hg19uniqmap.zip"
    conda:
        "../envs/basics.yaml"
    params:
        link=config["ref"]["hg19"]["uniqmapLink"]
    shell:
        "curl -fsSL {params.link} > {output}"

rule unzip_hg19_uniqmap:
    input:
        "data/2019-07-26/hg19uniqmap.zip"
    output:
        directory("data/2019-07-26/hg19-50nt-uniqmap")
    log:
        "logs/data/hg19_uniqmap.log"
    conda:
        "../envs/basics.yaml"
    shell:
        "unzip {input} -d data/2019-07-26"
