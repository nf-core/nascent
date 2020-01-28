rule bowtie2:
    """
    Aligns fastq samples with bowtie2 using reference genome
    """
    input:
        sample=get_fastqs,
        ref="data/2018-06-24/{genome}/genome.fa",
    output:
        "results/2018-10-04/{genome}/{sample}/Aligned.out.bam"
    log:
        "logs/{genome}/bowtie2/{sample}.log"
    params:
        # index=config["ref"]["index"],
        index= lambda w: config["ref"]["{}".format(w.genome)]["index"],
        extra=config["bowtie2"]["extra"]
    threads: 32
    wrapper:
        "0.38.0/bio/bowtie2/align"

# FIXME
# rule star_index:
#     """
#     Creates a star index
#     """
#     input:
#         fasta = "data/2018-10-16/fasta/{genome}/genome.fa",
#         gtf = "data/2018-11-09/{genome}/genes.gtf"
#     output:
#         directory("data/2018-10-16/{genome}/star")
#     message:
#         "Testing STAR index"
#     threads: 32
#     params:
#         extra = config["STAR"]["index"]["extra"]
#     log:
#         "logs/star_index_{genome}.log"
#     wrapper:
#         "0.40.0/bio/star/index"

# rule star_se:
#     """
#     Aligns fastq files with star
#     Requires more than 31G of memory
#     """
#     input:
#         fq1=get_fastqs,
#         starIndex = "data/2018-10-16/{genome}/star"
#     output:
#         "results/2018-10-04/{genome}/{sample}/Aligned.out.bam"
#     log:
#         "logs/{genome}/star/{sample}.log"
#     params:
#         # path to STAR reference genome index
#         index= lambda w: config["ref"]["{}".format(w.genome)]["starIndex"],
#         # optional parameters
#         extra = config["STAR"]["align"]["extra"]
#     threads: 32
#     resources:
#         mem_mb=32000
#     wrapper:
#         "0.40.0/bio/star/align"
