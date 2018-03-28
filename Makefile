# Makefile for IMR90 Analysis

# Run fastqc on all samples
getIMR90Data:

qualityCheck: getIMR90Data
	./results/2017-06-22/bin/runFastqc.sh

getRefGenome:

refGenomeAlignment: getRefGenome getIMR90Data
	./results/2017-06-27/bin/bowtie4all.sh
	./results/2017-06-27/bin/moveFastq.sh

createHomerTags: refGenomeAlignment

createHomerPeaks: createHomerTags

predicteRNAs: createHomerPeaks

getGMData:

compareIMR90toGM: predicteRNAs getGMData

covertSam:
	./results/2018-02-06/runall.sh

differentialAnalysis: compareIMR90toGM covertSam
# Run Create Coverage
	./results/2018-03-02/runall.sh

CompressAll:
	tar -zcvf /results/2017-06-27/shIMR90_sam.tar.gz /results/2017-06-27/sam/
