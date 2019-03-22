# nf-core/nascent: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [fastq-dump](#fastqdump) - if needed, extract the fastq file[s] from a sample
* [SeqKit/bbduk](#seqkitbbduk) - trim reads and remove adapters
* [FastQC](#fastqc) - read quality control
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline
* [HISAT2](#hisat2) - map reads to the reference genome
* [Samtools](#samtools) - convert the mapped reads as SAM files to BAM format
* [Preseq](#preseq) - estimate complexity of the sample
* [RSeQC](#rseqc) - analyze read distributions
* [Pileup](#pileup) - analyze coverage
* [bedtools](#bedtools) - create both normalized and non-normalized coverage files in bedGraph format
* [igvtools](#igvtools) - create compressed files to visualize the sample in the Integrative Genomics Viewer ([IGV](http://software.broadinstitute.org/software/igv/home))


## fastqdump
[fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) decompresses an SRR file obtained from the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) database. This will produce one or two fastq files (in the case of paired-end reads).

**Output directory: `results/fastq-dump`**

* `sample.fastq`
  * FastQ file to process, from the corresponding sample.


## seqkitbbduk
[SeqKit](https://bioinf.shenwei.me/seqkit/) is a toolkit for fasta and fastq file manipulation, used in the pipeline for xxxxxxxxxxxxxxxxxxxxxxxxx. [BBDuk](https://www.geneious.com/plugins/bbduk/) is an adapter trimming tool used to leave only the useful part of a sequenced read.

**Output directory: `results/bbduk`**

* `sample.trim.fastq`
  * Trimmed FastQ file for each sample.


## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/qc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info


## hisat2
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) is a sequence alignment tool to map the trimmed sequenced reads to the corresponding reference genome. Due to their size, the resulting sam files are not conserved after the pipeline has completed execution.

If the necessary indices for mapping are not provided/present, a separate process will build them first. This step can take a few minutes, however it should only be executed once.

**Output directory: none**



## samtools
[Samtools](http://www.htslib.org/) is a suite of tools to handle format conversions, among other things, for high-throughput sequencing data. We also use Samtools to generate the list of chromosome sizes, if not provided for the desired reference genome.

**Output directory: `results/mapped/bams`**

* `sample.trim.sorted.bam`
  * Mapped sample in BAM format
* `sample.trim.sorted.bam.bai`
  * Index for the `sample.trim.sorted.bam` mapped sample in BAM format

**Output directory: `results/qc/mapstats`**

* `sample.trim.sorted.bam.flagstat`
  * xxxxxxxxxxxxxxxxxxx
* `sample.trim.sorted.bam.millionsmapped`
  * xxxxxxxxxxxxxxxxxxx


## preseq
[Preseq](http://smithlabresearch.org/software/preseq/) plots the estimated complexity of a sample, and estimates future yields for complexity if the sample is sequenced at higher read depths.

**Output directory: `results/qc/preseq`**

* `sample.trim.c_curve.txt`
  * xxxxxxxxxxxxxxxxx
* `sample.trim.lc_extrap.txt`
  * xxxxxxxxxxxxxxxxx


## rseqc
[RSeQC](http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/) provides a number of useful modules that can comprehensively evaluate high throughput sequence data. We use it on this pipeline to analyze read distributions.

**Output directory: `results/qc/rseqc`**

* `sample.trim.read_dist.txt`
  * xxxxxxxxxxxxxxxxx


## pileup
[Pileup](xxxxxxxxxxxxxxxx) analyzes the sequencing coverage for each sample.

**Output directory: `results/qc/pileup`**

* `sample.trim.coverage.hist.txt`
  * xxxxxxxxxxxxxxxxx
* `sample.trim.coverage.stats.txt`
  * xxxxxxxxxxxxxxxxx


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

