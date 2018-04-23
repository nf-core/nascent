```
IMR90.Snakefile
Edmund Miller

Predict eRNA Transcripts from GRO-Seq data
``` 
from os.path import join, basename, dirname 

IMR0h IMR30min IMR1h IMR2h IMR4h IMR6h IMR12h IMR24h   

rule all:

# rule getIMR90Data:
#     output:
#         "data/2017-06-21/{IMR90}.fastq"
#     shell:    
#         # TODO
#         "wget https://functionalgenomics.org/data/gro-seq/IMR90/{IMR90}.fastq"

rule qualityCheck: 
    input: 
        "data/2017-06-21/{IMR90}.fastq"
    output:
        "results/2017-06-22/{IMR90}_fastqc.zip"
    shell:
        "fastqc {input} -o={output} -t=6"
    
rule getRefGenome:
    output:
        "data/2017-06-27/hg19.zip"
    shell:
        "wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip --output-document={output}"
        
rule buildRefGenome:
    input:
        "data/2017-06-27/hg19.zip"
    output:
        "hg19.1.bt2",       
        "hg19.2.bt2",       
        "hg19.3.bt2",       
        "hg19.4.bt2"
    threads: 4
    log:
        join(dirname(CDNA), 'kallisto/kallisto.index.log')
    run:
        shell('unzip {input}')
        shell('bash make_hg19.sh 2> {log}')
# rule refGenomeAlignment: getRefGenome getIMR90Data
# 	./results/2017-06-27/bin/bowtie4all.sh
# 	./results/2017-06-27/bin/moveFastq.sh

# rule createHomerTags: refGenomeAlignment

# rule createHomerPeaks: createHomerTags

# rule predicteRNAs: createHomerPeaks

# rule getGMData:

# rule compareIMR90toGM: predicteRNAs getGMData

# rule covertSam:
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
