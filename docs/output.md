# nf-core/nascent: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [FastQC](#fastqc) - Raw read QC
  - [fastp](#fastp) - Adapter and quality trimming
- [Alignment](#alignment)
  - [bwa](#bwa) - TODO
  - [bwa-mem2](#bwa-mem2) - TODO
  - [DRAGMAP](#dragmap) - TODO
- [Alignment post-processing](#alignment-post-processing)
  - [SAMtools](#samtools) - Sort and index alignments
  - [UMI-tools dedup](#umi-tools-dedup) - UMI-based deduplication
  - [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking
- [Quality control](#quality-control)
  - [RSeQC](#rseqc) - Various RNA-seq QC metrics
  - [Preseq](#preseq) - Estimation of library complexity
  - [BBMap](#bbmap) - TODO
- [Coverage Graphs](#coverage-graphs)
  - [BEDTools Genomecov](#bedtools-genomcov) - Create bigWig coverage files
  - [deepTools bamcoverage](#deeptools-bamcoverage) - TODO
- [Transcript Identification](#transcript-identification)
  - [GroHMM](#grohmm) - Predicts transcripts from aligned GROSeq data in the form of bed files.
  - [HOMER](#homer) - TODO
  - [PINTS](#pints) - TODO
  - [BEDTools Insersect](#bedtools-intersect) - TODO
- [Quantification](#quantification)
  - [featureCounts](#featurecounts) - Read counting relative to gene biotype
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `<flowcell_id>/`
  - `*.fastp.html`: Trimming report in html format.
  - `*.fastp.json`: Trimming report in json format.
  - `*.fastp.log`: Trimming log file.

</details>

[fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast, all-in-one preprocessing for FastQ files. It has been developed in C++ with multithreading support to achieve higher performance. fastp is used in this pipeline for standard adapter trimming and quality filtering.

![MultiQC - fastp filtered reads plot](images/mqc_fastp_plot.png)

## Alignment

### bwa

<!-- TODO Add files -->

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome. The aligned reads are then coordinate-sorted is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

### bwa-mem2

<!-- TODO Add files -->

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) is a software package for mapping low-divergent sequences against a large reference genome.The aligned reads are then coordinate-sorted with [samtools](https://www.htslib.org/doc/samtools.html).

### DRAGMAP

<!-- TODO Add files -->

[DragMap](https://github.com/Illumina/dragmap) is an open-source software implementation of the DRAGEN mapper, which the Illumina team created so that we would have an open-source way to produce the same results as their proprietary DRAGEN hardware. The aligned reads are then coordinate-sorted with [samtools](https://www.htslib.org/doc/samtools.html).

## Alignment post-processing

### SAMtools

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/`
  - `<SAMPLE>.sorted.bam`: If `--save_align_intermeds` is specified the original coordinate sorted BAM file containing read alignments will be placed in this directory.
  - `<SAMPLE>.sorted.bam.bai`: If `--save_align_intermeds` is specified the BAI index file for the original coordinate sorted BAM file will be placed in this directory.
  - `<SAMPLE>.sorted.bam.csi`: If `--save_align_intermeds --bam_csi_index` is specified the CSI index file for the original coordinate sorted BAM file will be placed in this directory.
- `<ALIGNER>/samtools_stats/`
  - SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_mapped.png)

![MultiQC - SAMtools mapped reads per contig plot](images/mqc_samtools_idxstats.png)

### UMI-tools dedup

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/`
  - `<SAMPLE>.umi_dedup.sorted.bam`: If `--save_umi_intermeds` is specified the UMI deduplicated, coordinate sorted BAM file containing read alignments will be placed in this directory.
  - `<SAMPLE>.umi_dedup.sorted.bam.bai`: If `--save_umi_intermeds` is specified the BAI index file for the UMI deduplicated, coordinate sorted BAM file will be placed in this directory.
  - `<SAMPLE>.umi_dedup.sorted.bam.csi`: If `--save_umi_intermeds --bam_csi_index` is specified the CSI index file for the UMI deduplicated, coordinate sorted BAM file will be placed in this directory.
- `<ALIGNER>/umitools/`
  - `*_edit_distance.tsv`: Reports the (binned) average edit distance between the UMIs at each position.
  - `*_per_umi.tsv`: UMI-level summary statistics.
  - `*_per_umi_per_position.tsv`: Tabulates the counts for unique combinations of UMI and position.

The content of the files above is explained in more detail in the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html#dedup-specific-options).

</details>

After extracting the UMI information from the read sequence (see [UMI-tools extract](#umi-tools-extract)), the second step in the removal of UMI barcodes involves deduplicating the reads based on both mapping and UMI barcode information using the UMI-tools `dedup` command. This will generate a filtered BAM file after the removal of PCR duplicates.

### picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/`
  - `<SAMPLE>.markdup.sorted.bam`: Coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM file and so will be saved by default in the results directory.
  - `<SAMPLE>.markdup.sorted.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory.
  - `<SAMPLE>.markdup.sorted.bam.csi`: CSI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory. Only generated if `--bam_csi_index` is specified as a parameter.
- `<ALIGNER>/samtools_stats/`
  - SAMtools `<SAMPLE>.markdup.sorted.bam.flagstat`, `<SAMPLE>.markdup.sorted.bam.idxstats` and `<SAMPLE>.markdup.sorted.bam.stats` files generated from the duplicate marked alignment files.
- `<ALIGNER>/picard_metrics/`
  - `<SAMPLE>.markdup.sorted.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

Unless you are using [UMIs](https://emea.illumina.com/science/sequencing-method-explorer/kits-and-arrays/umi.html) it is not possible to establish whether the fragments you have sequenced from your sample were derived via true biological duplication (i.e. sequencing independent template fragments) or as a result of PCR biases introduced during the library preparation. By default, the pipeline uses [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) to _mark_ the duplicate reads identified amongst the alignments to allow you to guage the overall level of duplication in your samples. However, for RNA-seq data it is not recommended to physically remove duplicate reads from the alignments (unless you are using UMIs) because you expect a significant level of true biological duplication that arises from the same fragments being sequenced from for example highly expressed genes. You can skip this step via the `--skip_markduplicates` parameter.

![MultiQC - Picard MarkDuplicates metrics plot](images/mqc_picard_markduplicates.png)

## Quality control

### RSeQC

[RSeQC](<(http://rseqc.sourceforge.net/)>) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. You can tweak the supported scripts you would like to run by adjusting the `--rseqc_modules` parameter which by default will run all of the following: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

### Preseq

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/preseq/`
  - `*.lc_extrap.txt`: Preseq expected future yield file.
- `<ALIGNER>/preseq/log/`
  - `*.command.log`: Standard error output from command.

</details>

The [Preseq](http://smithlabresearch.org/software/preseq/) package is aimed at predicting and estimating the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. The estimates can then be used to examine the utility of further sequencing, optimize the sequencing depth, or to screen multiple libraries to avoid low complexity samples. A shallow curve indicates that the library has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads. Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

![MultiQC - Preseq library complexity plot](images/mqc_preseq_plot.png)

### BBMap

<!-- TODO -->

## Coverage Graphs

### BEDTools Genomecov

<!-- TODO -->

### deepTools bamcoverage

<!-- TODO -->

## Transcript Identification

### HOMER

<!-- TODO -->

### PINTS

<!-- TODO -->

### GroHMM

<details markdown="1">
<summary>Output files</summary>

- `grohmm/`
  - `*.eval.txt`: Evaluation of HMM Annotations
  - `*.final.transcripts.bed`: Predicted transcripts
  - `*.tdFinal.txt`: Final quality metrics
  - `*.tdplot_mqc.jpg`: TD plot included in multiqc
  - `*.transcripts.txt`: Predicted transcripts in txt form
  - `*.tuning.csv`: The tuning csv that was used

</details>

[GroHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html) is a computational tool for identifying unannotated and cell type-specific transcription units from GRO-seq data. The pipeline will predict, and then repair transcripts based on known errors to generate a final set of transcripts (in the form of a bed file) for further analysis.
By default, tuning will be performed by inputting a preset comma-separated values file with two columns, each identifying tuning parameters - LtProbB and UTS. These refer to the log-transformed transition probability of switching from transcribed state to non-transcribed state and variance of the emission probability for reads in the non-transcribed state, respectively. The output of the tuning file, also a comma-separated values file, will list out the sum of errors and error rate per called transcript, which will enable Nextflow to specify optimal UTS and LtProbB values for the subsequent transcript identification step. The user may also choose to provide their own list of hold-out parameters to test (in the format of a .csv file), or skip the tuning process altogether due to time constraints. If the tuning process is skipped ('--skip_tuning') then the user may indicate the specific holdout parameters to use ('--uts' and '--ltprobb') or choose to use the default parameters.
The transcript calling step will use the two-state hidden Markov model (HMM) which GroHMM employs in order to identify boundaries of transcription across the genome in a de-novo manner. The output is a .bed file of transcripts used in downstream analysis.

For more information about how to use GROHMM, see the [tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/groHMM/inst/doc/groHMM.pdf) or [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/groHMM/man/groHMM.pdf).

## Quatification

### featureCounts

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/featurecounts/`
  - `*.featureCounts.txt`: featureCounts biotype-level quantification results for each sample.
  - `*.featureCounts.txt.summary`: featureCounts summary file containing overall statistics about the counts.
  - `*_mqc.tsv`: MultiQC custom content files used to plot biotypes in report.

</details>

[featureCounts](http://bioinf.wehi.edu.au/featureCounts/) from the [Subread](http://subread.sourceforge.net/) package is a quantification tool used to summarise the mapped read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations. We can also use featureCounts to count overlaps with different classes of genomic features. This provides an additional QC to check which features are most abundant in the sample, and to highlight potential problems such as rRNA contamination.

<!-- TODO Add example plot -->

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
