# Makefile for IMR90 Analysis

# Run fastqc on all samples
getIMR90Data:

qualityCheck: getIMR90Data

getRefGenome:

refGenomeAlignment: getRefGenome getIMR90Data

createHomerTags: refGenomeAlignment

createHomerPeaks: createHomerTags

predicteRNAs: createHomerPeaks

getGMData:

compareIMR90toGM: predicteRNAs getGMData

differentialAnalysis: compareIMR90toGM
	# Run Create Coverage
	results/2018-03-02/runall.sh
