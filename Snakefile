'''
IMR90.Snakefile
Edmund Miller

Predict eRNA Transcripts from GRO-Seq data
'''
from os.path import join, basename, dirname 

IMR90 = ["IMR0h", "IMR30min", "IMR1h", "IMR2h", "IMR4h", "IMR6h", "IMR12h", "IMR24h"]   
configfile: 'config.yml'

# hg19 = config['hg19']

rule all:
    input:
        "results/2017-08-02/GroseqIMR90peak.gtf",
        "results/2018-01-25/eRNA_IMR90_GM_hg19_Overlap.bed"
# rule report:
#     input:
#         "results/2017-07-27/All_together/IMR90tsspeak.txt"
#     output:
#         "report.html"
#     run:
#         from snakemake.utils import report
#         with open(input[0]) as vcf:
#             n_calls = sum(1 for l in vcf if not l.startswith("#"))

#         report("""
#         An example variant calling workflow
#         ===================================

#         Reads were mapped to the Yeast
#         reference genome and variants were called jointly with
#         SAMtools/BCFtools.

#         This resulted in {n_calls} variants (see Table T1_).
#         """, output[0], T1=input[0])
rule IMR90_Data:
    output:
        "data/2017-06-21/{IMR90}.fastq"
    run:    
        # TODO https
        shell('wget http://functionalgenomics.org/data/gro-seq/IMR90_GROseq.tar.gz')
        shell('tar -xvf IMR90_GROseq.tar.gz data/2017-06-21/.')

rule fastqc: 
    input: 
        "data/2017-06-21/{IMR90}.fastq"
    output:
        "results/2017-06-22/{IMR90}_fastqc.zip"
    shell:
        "fastqc {input} -o={output} -t=6"
    
rule Bowtie2_Reference_Genome:
    output:
        "data/2017-06-27/hg19.zip"
    shell:
        "wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip --output-document={output}"
        
rule Build_Ref_Genome:
    input:
        "data/2017-06-27/hg19.zip"
    output:
        "hg19.1.bt2",       
        "hg19.2.bt2",       
        "hg19.3.bt2",       
        "hg19.4.bt2"
    threads: 4
    log:
        #join(dirname(hg19), 'bt2.log')
        "data/2017-06-27/buildRefGenome.log"
    run:
        shell('unzip {input}')
        shell('bash make_hg19.sh 2> {log}')

# FIXME Hard coded the bowtie2 index
rule Bowtie2:
    input:
        "data/2017-06-21/{IMR90}.fastq",
        index = "data/2017-06-27/hg19.1.bt2"
    output:
        "results/2017-06-27/{IMR90}.sam"
    log:
        "results/2017-06-27/bowtie2Alignment.log"
    # FIXME
    # conda:
    #     "envs/bowtie2.yaml"
    threads:4
    run: 
        shell('bowtie2'
              ' --time'
              ' --threads={threads}'
              ' --index=data/2017-06-27/hg19.1.bt2' #{input.index}'
              ' -- very-sensitive'
              '-x {IMR90}'
              '-S {output}')

rule makeTagDirectory:
    input:
        "results/2017-06-27/{IMR90}.sam"
    output:
        "results/2017-07-27/Tags/{IMR90}/"
    shell:
        "makeTagDirectory tags/{input} -genome hg19 -checkGC {output}"
        
rule AllTogether_makeTagDirectory:
    input:
        "results/2017-06-27/*.sam"
    output:
        "results/2017-07-27/All_together/"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"

rule findPeaks:
    input:
        "results/2017-07-27/All_together/"
    output:
        "results/2017-08-02/GroseqIMR90peak.gtf"
    shell:
        "findPeaks {input} -o {output} -style groseq"

# FIXME
rule pos2bed: 
    input:
        "results/2017-08-02/GroseqIMR90peak.gtf"
    output:
        "results/2017-08-02/GroseqIMR90peak.bed"
    run:
     shell('pos2bed.pl {input} > {output}')
     # shell('bedtools intersect -a GroseqIMR90peak.bed -b slop-refseqhg19.bed -v > GroseqIMR90nogenes.bed')

# FIXME
rule getRefSeqhg19:     
    output: 
        "data/2017-08-02/refseqhg19.bed"

# FIXME
rule slopRefSeq:
    input:
        "data/2017-08-02/refseqhg19.bed"
    output:
        "data/2017-08-02/sloprefseqhg19.bed"

# FIXME
rule RemoveGenes:
    input:
        IMR = "results/2017-08-02/GroseqIMR90peak.bed",
        refseq = "data/2017-08-02/sloprefseqhg19.bed"
    output:
        "results/2017-08-02/GroseqIMR90nogenes.bed"
    run:
     shell('bedtools intersect -a {input.IMR} -b {refseq} -v > {output}')

# FIXME     
rule getHistones:
    output:
        "results/2017-10-23/E017-H3K27ac.pval.signal.bigwig",
        "results/2017-10-23/E017-H3K4me1.pval.signal.bigwig"

# FIXME
rule Histones:
    input:
        IMR90 = "results/2017-08-02/GroseqIMR90nogenes.bed",
        H3K27ac = "results/2017-10-23/E017-H3K27ac.pval.signal.bigwig",
        H3K4me1 = "results/2017-10-23/E017-H3K4me1.pval.signal.bigwig"
    output:
        "results/2017-11-01/eRNA_IMR90_hg19.bed"

# FIXME
rule GMData:
    output:
        "data/2018-01-25/eRNA_GM_hg19.bed"
    shell:
        "http://functionalgenomics.org/data/gro-seq/GM/eRNA_GM_hg19.bed {output}"

# FIXME
rule intersectBed:
    input:
        eGM = "data/2018-01-25/eRNA_GM_hg19.bed",
        eIMR = "results/2017-11-01/eRNA_IMR90_hg19.bed"
    output:
        "results/2018-01-25/eRNA_IMR90_GM_hg19_Overlap.bed"
    run:
        shell('intersectBed -wo -a {input.eIMR} -b {input.eGM} > {output}')

# Sam:
# 	./results/2018-02-06/bin/covert2bam.sh
# 	./results/2018-02-06/bin/covert2bed.sh

# rule differentialAnalysis: compareIMR90toGM covertSam
# # Run Create Coverage
# 	./results/2018-03-02/bin/run_Unique.sh
# 	./results/2018-03-02/bin/run_Unique_Overlap.sh
# 	./results/2018-03-02/bin/run_Overlap.sh
# 	./results/2018-03-02/bin/run_GM_Overlap.sh
# 	./results/2018-03-02/bin/run_GM.sh

# rule CompressAll:
# 	tar -zcvf /results/2017-06-27/shIMR90_sam.tar.gz /results/2017-06-27/sam/
