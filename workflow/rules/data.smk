# TODO move path to config
rule AWS_iGenomes:
    output:
        "scripts/aws-igenomes.sh"
    priority: 50
    shell:
       "curl -fsSL https://ewels.github.io/AWS-iGenomes/aws-igenomes.sh > {output} && chmod +x {output}"

# Datasets TODO

# TODO move path to config
rule GM_download:
    output:
        expand("data/2018-06-23/{unit}.fastq",unit=GM_SAMPLES),

# TODO move path to config
rule IMR_download:
    output:
        expand("data/2018-06-23/{unit}.fastq",unit=IMR_SAMPLES),


#############################
# BowTie2 Reference Genomes #
#############################

# TODO move path to config
rule bowtie2_reference_Genome:
    output:
        "data/2018-06-24/{unit}/genome.fa",
    params:
        script="scripts/aws-igenomes.sh",
        genome="Homo_sapiens",
        source="UCSC",
        build="{unit}",
        typeOf="bowtie2",
        outDir="data/2018-06-24/{unit}",
    conda:
        "../envs/awscli.yaml"
    priority: 50
    shell:
        "{params.script} -g {params.genome} -s {params.source} "
        "-b {params.build} -t {params.typeOf} -o {params.outDir}"

rule genome_fasta:
    output:
        "data/2018-10-16/fasta/{unit}/genome.fa",
    params:
        script="scripts/aws-igenomes.sh",
        genome="Homo_sapiens",
        source="UCSC",
        build="{unit}",
        typeOf="fasta",
        outDir="data/2018-10-16/fasta/{unit}",
    conda:
        "../envs/awscli.yaml"
    priority: 50
    log:
        "logs/{unit}/genome.log"
    shell:
        "{params.script} -g {params.genome} -s {params.source} "
        "-b {params.build} -t {params.typeOf} -o {params.outDir}"

# TODO move path to config
rule fai:
    input:
        "data/2018-06-24/{unit}/genome.fa",
    output:
        "data/2018-06-24/{unit}/genome.fa.fai",
    log:
        "logs/{unit}_fai.log"
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} > {output}"

# TODO move path to config
rule chrom_len:
    input:
        "data/2018-06-24/{unit}/genome.fa.fai",
    output:
        "data/2018-06-24/{unit}/chrom.sizes",
    log:
        "logs/{unit}_chrom.sizes.log"
    threads: 4
    # conda:
    #     "../envs/samtools.yaml"
    shell:
        "cut -f 1,2 {input} > {output}"

# TODO move path to config
rule hg18_UCSC_chrom_size:
    output:
        "data/2019-09-25/hg18.chrom.sizes"
    shell:
       "curl -fsSL https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > {output}"

# TODO move path to config
rule hg19_UCSC_chrom_size:
    output:
        "data/2019-09-25/hg19.chrom.sizes"
    shell:
       "curl -fsSL https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > {output}"

##########
# RefSeq #
##########
# TODO move path to config

rule refSeq_gtf:
    output:
        "data/2018-11-09/{unit}/genes.gtf",
    params:
        script="scripts/aws-igenomes.sh",
        genome="Homo_sapiens",
        source="UCSC",
        build="{unit}",
        typeOf="gtf",
        outDir="data/2018-11-09/{unit}/",
    conda:
        "../envs/awscli.yaml"
    shell:
        "{params.script} -g {params.genome} -s {params.source} "
        "-b {params.build} -t {params.typeOf} -o {params.outDir}"

rule refSeq_bed:
    output:
        "data/2018-11-09/{unit}/genes.bed",
    params:
        script="scripts/aws-igenomes.sh",
        genome="Homo_sapiens",
        source="UCSC",
        build="{unit}",
        typeOf="bed12",
        outDir="data/2018-11-09/{unit}/",
    conda:
        "../envs/awscli.yaml"
    shell:
        "{params.script} -g {params.genome} -s {params.source} "
        "-b {params.build} -t {params.typeOf} -o {params.outDir}"

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
