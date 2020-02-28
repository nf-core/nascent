# nf-core/nascent: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [fastq-dump](#fastqdump) - if needed, extract the fastq file[s] from a sample
* [SeqKit/bbduk](#seqkitbbduk) - flip reads (experiment specific) & trim reads for adapters/quality/length
* [FastQC](#fastqc) - read quality control
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [HISAT2](#hisat2) - map reads to the reference genome
* [Samtools](#samtools) - convert the mapped reads as SAM files to BAM format
* [Preseq](#preseq) - estimate complexity of the sample
* [RSeQC](#rseqc) - analyze read distributions, infer experiment (SE/PE, whether reads need to be flipped), & read duplication
* [BBMap](#pileup) - analyze coverage
* [bedtools](#bedtools) - create both normalized and non-normalized coverage files in bedGraph format
* [igvtools](#igvtools) - create compressed files to visualize the sample in the Integrative Genomics Viewer ([IGV](http://software.broadinstitute.org/software/igv/home))


## fastqdump
[fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) decompresses an SRR file obtained from the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) database. This will produce one or two fastq files (in the case of paired-end reads).

**Output directory: `results/fastq-dump`**

* `sample.fastq`
  * FastQ file to process, from the corresponding sample.


## seqkit & bbduk
[SeqKit](https://bioinf.shenwei.me/seqkit/) is a toolkit for fasta and fastq file manipulation, used in the pipeline if the positive/negative strands need to be flipped (dependent on library prep protocol). [BBDuk](https://www.geneious.com/plugins/bbduk/) is trimming tool used to filter reads for adapters, read quality, and overall length after adapter removal.

**Output directory: `results/bbduk, qc/trimstats`**

* `sample.trim.fastq`
  * Trimmed FastQ file for each sample.
* `{refstats,trimstats,ehist}.txt`
  * Trimming details including adapters removed, percentages of reads removed that did not meet minimum quality/length


## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows both untrimmed and trimmed reads.

**Output directory: `results/qc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files & trimmed fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## MultiQC
[MultiQC](https://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info)


## hisat2
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) is a sequence alignment tool to map the trimmed sequenced reads to the corresponding reference genome. Due to their size, the resulting sam files are not conserved after the pipeline has completed execution.

If the necessary indices for mapping are not provided/present, a separate process will build them first. This step can take a few minutes, however it should only be executed once.

## samtools
[Samtools](http://www.htslib.org/) is a suite of tools to handle format conversions, among other things, for high-throughput sequencing data. We also use Samtools to generate the list of chromosome sizes, if not provided for the desired reference genome.

**Output directory: `results/mapped/bams`**

* `sample.trim.sorted.bam`
  * Mapped sample in BAM format
* `sample.trim.sorted.bam.bai`
  * Index for the `sample.trim.sorted.bam` mapped sample in BAM format

**Output directory: `results/qc/mapstats`**

* `sample.trim.sorted.bam.flagstat`
  * Overall mapping statistics
* `sample.trim.sorted.bam.millionsmapped`
  * File that contains number of uniquely mapped reads (not total multi-mapped). Used in normalization


## preseq
[Preseq](http://smithlabresearch.org/software/preseq/) plots the estimated complexity of a sample, and estimates future yields for complexity if the sample is sequenced at higher read depths.

**Output directory: `results/qc/preseq`**

* `sample.trim.c_curve.txt`
  * Curve generated based on number of unique reads vs. total reads sequenced
* `sample.trim.lc_extrap.txt`
  * Extrapolation of the c_curve that attempts to model the predicted number of unique reads if the sample was seqeunced to a greater depth


## rseqc
[RSeQC](http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/) provides a number of useful modules that can comprehensively evaluate high throughput sequence data. We use it on this pipeline to analyze read distributions.

**Output directory: `results/qc/rseqc`**

* `sample.trim.read_dist.txt`
  * Relative distribution of reads relative to a gene reference file


## pileup
[BBMap](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh) includes a tool called `pileup`, which analyzes the sequencing coverage for each sample.

**Output directory: `results/qc/pileup`**

* `sample.trim.coverage.hist.txt`
  * Histogram of read coverage over each chromosome
* `sample.trim.coverage.stats.txt`
  * Coverage stats broken down by chromosome including %GC, pos/neg read coverage, total coverage, etc.


## bedtools
[bedtools](https://bedtools.readthedocs.io/en/latest/) is an extensive toolkit for BED and bedGraph format manipulation, like sorting, intersecting and joining these files. The files produced here are useful to be processed later using [Tfit](https://github.com/Dowell-Lab/Tfit) or [dReg](https://github.com/Danko-Lab/dREG) to find regions of active transcription, and transcription regulatory elements.

**Output directory: `results/mapped/bedgraphs`**

* `sample.trim.bedGraph`
  * Sample coverage file in bedGraph format
* `sample.trim.pos.bedGraph`
  * Sample coverage file (positive strand only) in bedGraph format
* `sample.trim.neg.bedGraph`
  * Sample coverage file (negative strand only) in bedGraph format

**Output directory: `results/mapped/rcc_bedgraphs`**

* `sample.trim.rcc.bedGraph`
  * Normalized sample coverage file in bedGraph format
* `sample.pos.trim.rcc.bedGraph`
  * Normalized sample coverage file (positive strand only) in bedGraph format
* `sample.neg.trim.rcc.bedGraph`
  * Normalized sample coverage file (negative strand only) in bedGraph format

**Output directory: `results/mapped/dreg_input`**

* `sample.trim.pos.rcc.bw`
  * Sample coverage file (positive strand only) in BigWig format
* `sample.trim.neg.rcc.bw`
  * Sample coverage file (negative strand only) in BigWig format

**Output directory: `results/mapped/rcc_bigwig`**

* `sample.trim.pos.rcc.bw`
  * Normalized sample coverage file (positive strand only) in BigWig format
* `sample.trim.neg.rcc.bw`
  * Normalized sample coverage file (negative strand only) in BigWig format


## igvtools
[igvtools](https://software.broadinstitute.org/software/igv/igvtools) is a commandline tool we use to produce a compressed version of the sample coverage file in order to visualize it on IGV more efficiently (with a significantly smaller memory footprint).

**Output directory: `results/mapped/tdfs`**

* `sample.trim.rpkm.tdf`
  * Sample coverage file in TDF format

