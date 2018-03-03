#https://stackoverflow.com/questions/26424594/how-to-get-a-makefile-to-run-all-commands-regardless-of-targets-or-dependencies
objects = IMR0H_eRNA_IMR90_unique_coverage.bam \
IMR30min_eRNA_IMR90_unique_coverage.bam \
IMR1h_eRNA_IMR90_unique_coverage.bam \
IMR2h_eRNA_IMR90_unique_coverage.bam \
IMR4h_eRNA_IMR90_unique_coverage.bam \
IMR6h_eRNA_IMR90_unique_coverage.bam \
IMR12h_eRNA_IMR90_unique_coverage.bam\
IMR24h_eRNA_IMR90_unique_coverage.bam

all: $(objects)

%_eRNA_IMR90_unique_coverage.bam: ../2018-02-06/bam/%.bam ../2018-01-29/eRNA_IMR90_hg19_No_Overlap.bed
	# Run coverageBed on all available bam files
	coverageBed -d -header -sorted -a *h.bam -b eRNA > $@
