<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nascent_logo_dark.png">
    <img alt="nf-core/nascent" src="docs/images/nf-core-nascent_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/nascent/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/nascent/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nascent/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nascent/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nascent/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7245273-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7245273)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/nascent)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nascent-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nascent)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/nascent** is a bioinformatics best-practice analysis pipeline for nascent transcript (NT) and Transcriptional Start Site (TSS) assays.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/nascent/results).

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))
3. Alignment
   1. [`bwa`](https://bio-bwa.sourceforge.net/)
   2. [`bwamem2`](https://github.com/bwa-mem2/bwa-mem2)
   3. [`DRAGMAP`](https://github.com/Illumina/DRAGMAP)
4. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
6. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
7. Quality Control
   1. [`RSeQC`](https://rseqc.sourceforge.net/index.html) - Various RNA-seq QC metrics
   2. [`Preseq`](http://smithlabresearch.org/software/preseq/) - Estimation of library complexity
   3. [`BBMap`](https://sourceforge.net/projects/bbmap/) - Analyzes the sequencing coverage
8. Coverage Graphs
   1. Create bedGraph coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/)
   2. Create bigWig coverage files ([`deeptools`](https://deeptools.readthedocs.io/en/develop/))
9. Transcript identification
   1. [`HOMER`](http://homer.ucsd.edu/)
   2. [`GroHMM`](https://bioconductor.org/packages/release/bioc/html/groHMM.html)
   3. [`PINTS`](https://pints.yulab.org/)
10. Quantification of Genes and Nascent Transcripts ([`featureCounts`](https://subread.sourceforge.net/featureCounts.html))
11. Aggregate report describing results and QC from the whole pipeline ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:


```csv title="samplesheet.csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

```bash
nextflow run nf-core/nascent \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/nascent/usage) and the [parameter documentation](https://nf-co.re/nascent/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/nascent/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/nascent/output).

## Credits

nf-core/nascent was originally written by Ignacio Tripodi ([@ignaciot](https://github.com/ignaciot)) and Margaret Gruca ([@magruca](https://github.com/magruca)).

The pipeline was re-written in Nextflow DSL2 by Edmund Miller ([@edmundmiller](https://github.com/edmundmiller)) and Sruthi Suresh ([@sruthipsuresh](https://github.com/sruthipsuresh)) from [The Functional Genomics Laboratory](https://taehoonkim.org/) at [The Univeristy of Texas at Dallas](https://www.utdallas.edu/)

We thank the following people for their extensive assistance in the development of this pipeline:

- [@apeltzer](https://github.com/apeltzer)
- [@ewels](https://github.com/ewels)
- [@drpatelh](https://github.com/drpatelh)
- [@pditommaso](https://github.com/pditommaso)
- [@FriederikeHanssen](https://github.com/FriederikeHanssen)
- [Tae Hoon Kim](https://github.com/taehoonkim-phd)
- [@easterwoods](https://github.com/easterwoods)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nascent` channel](https://nfcore.slack.com/channels/nascent) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/nascent for your analysis, please cite it using the following doi: [10.5281/zenodo.7245273](https://doi.org/10.5281/zenodo.7245273)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
