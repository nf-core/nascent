# nf-core/nascent: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [FastQC](#fastqc) - Raw read QC
  - [fastp](#fastp) - TODO
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

<!-- TODO -->

## Alignment

### bwa

<!-- TODO -->

### bwa-mem2

<!-- TODO -->

### DRAGMAP

<!-- TODO -->

## Alignment post-processing

### SAMtools

<!-- TODO -->

### UMI-tools dedup

<!-- TODO -->

### picard MarkDuplicates

<!-- TODO -->

## Quality control

### RSeQC

<!-- TODO -->

### Preseq

<!-- TODO -->

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

#### TODO: Add output files once full pipeline is run

</details>

[GroHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html) is a computational tool for identifying unannotated and cell type-specific transcription units from GRO-seq data. The pipeline will predict, and then repair transcripts based on known errors to generate a final set of transcripts (in the form of a bed file) for further analysis.
By default, tuning will be performed by inputting a preset comma-separated values file with two columns, each identifying tuning parameters - LtProbB and UTS. These refer to the log-transformed transition probability of switching from transcribed state to non-transcribed state and variance of the emission probability for reads in the non-transcribed state, respectively. The output of the tuning file, also a comma-separated values file, will list out the sum of errors and error rate per called transcript, which will enable Nextflow to specify optimal UTS and LtProbB values for the subsequent transcript identification step. The user may also choose to provide their own list of hold-out parameters to test (in the format of a .csv file), or skip the tuning process altogether due to time constraints. If the tuning process is skipped ('--skip_tuning') then the user may indicate the specific holdout parameters to use ('--uts' and '--ltprobb') or choose to use the default parameters.
The transcript calling step will use the two-state hidden Markov model (HMM) which GroHMM employs in order to identify boundaries of transcription across the genome in a de-novo manner. The output is a .bed file of transcripts used in downstream analysis.

For more information about how to use GROHMM, see the [tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/groHMM/inst/doc/groHMM.pdf) or [documentation](https://www.bioconductor.org/packages/release/bioc/manuals/groHMM/man/groHMM.pdf).

## Quatification

### featureCounts

<!-- TODO -->

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
