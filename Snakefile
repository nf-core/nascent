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

# rule IMR90_Data:
#     output:
#         "data/2017-06-21/{IMR90}.fastq"
#     run:    
#         # FIXME https
#         shell('wget http://functionalgenomics.org/data/gro-seq/IMR90_GROseq.tar.gz')
#         shell('tar -xvf IMR90_GROseq.tar.gz data/2017-06-21/.')

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
    run:
        shell('wget' 
            ' ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip' 
            ' --output-document={output}')
rule Build_bowtieRef_Genome:
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

rule Bowtie2:
    input:
        "data/2017-06-21/{IMR90}.fastq",
        index = "data/2017-06-27/hg19.1.bt2"
    output:
        "results/2017-06-27/{IMR90}.sam"
    log:
        "results/2017-06-27/bowtie2Alignment.log"
    # FIXME Hard coded the bowtie2 index
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

rule pos2bed: 
    input:
        "results/2017-08-02/GroseqIMR90peak.gtf"
    output:
        "results/2017-08-02/GroseqIMR90peak.bed"
    run:
        shell('pos2bed.pl {input} > {output}')

# FIXME
# Can't find where to redownload
rule RefSeqhg19:     
    output: 
        "data/2017-08-02/refseqhg19.bed"

# Not necessary just git fa file of lengths
# rule Genomehg19:
#     output:
#         tar = "data/2017-07-27/chromFa.tar.gz",
#     run:
#         shell('wget --timestamping'
#             ' ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
#             ' -O {output.tar}')
#         shell('tar xvzf {output.tar}')

rule slopRefSeq:
    input:
        refSeq="data/2017-08-02/refseqhg19.bed",
        chromLen="data/2017-07-27/genome.fa.fai"
    output:
        "results/2017-08-02/sloprefseqhg19.bed"
    run:
         shell('slopBed -i {input.refSeq}'
         ' -g {input.chromLen}'
         ' -l 1000 -r 10000 > {output}')

rule RemoveGenes:
    input:
        IMR = "results/2017-08-02/GroseqIMR90peak.bed",
        refseq = "results/2017-08-02/sloprefseqhg19.bed"
    output:
        "results/2017-08-02/GroseqIMR90nogenes.bed"
    log:
        "log/RemoveGenes.log"
    run:
     shell('bedtools intersect -a {input.IMR} -b {input.refseq} -v'
           ' | sort -k1,1 -k2,2n - > {output} 2> {log}')

# TODO pysam the unzipping and add conda support
rule wgetH3K27ac:
    output:
        "data/2017-10-23/GSM469967_UCSD.IMR90.H3K27ac.LL235.bed",
    log:
        "log/wgetH3K27ac.log"
    run:
        shell('wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM469nnn/GSM469967/suppl/GSM469967_UCSD.IMR90.H3K27ac.LL235.bed.gz -O data/2017-10-23/GSM469967_UCSD.IMR90.H3K27ac.LL235.bed.gz 2> {log}')
        shell('bgzip -d data/2017-10-23/GSM469967_UCSD.IMR90.H3K27ac.LL235.bed.gz'
           ' | sort -k1,1 -k2,2n -k3,3n -k4,4n - 2> {log}')

rule wgetH3K4me3:
    output:
        "data/2017-10-23/GSM469970_UCSD.IMR90.H3K4me3.LL221.bed"
    log:
        "log/wgetH3K4me3.log"
    run:
        shell('wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM469nnn/GSM469970/suppl/GSM469970_UCSD.IMR90.H3K4me3.LL221.bed.gz -O data/2017-10-23/GSM469970_UCSD.IMR90.H3K4me3.LL221.bed.gz 2> {log}')
        shell('bgzip -d data/2017-10-23/GSM469970_UCSD.IMR90.H3K4me3.LL221.bed.gz'
            ' | sort -k1,1 -k2,2n -k3,3n -k4,4n - 2> {log}')

# FIXME
rule HistonesIntersect:
    input:
        IMR90nogenes = "results/2017-08-02/GroseqIMR90nogenes.bed",
        H3K27ac = "data/2017-10-23/GSM469967_UCSD.IMR90.H3K27ac.LL235.bed",
        H3K4me3 = "data/2017-10-23/GSM469970_UCSD.IMR90.H3K4me3.LL221.bed"
    output:
        "results/2017-11-01/eRNA_IMR90_hg19.bed"
    log:
        "results/2017-11-01/HistonesIntersect.log"
    shell:
        "bedtools intersect -a {input.IMR90nogenes} -b {input.H3K27ac} {input.H3K4me3} -sorted -u > {output} 2> {log}"

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
