# ![nf-core/nascent](docs/images/nf-core-nascent_logo_light.png#gh-light-mode-only) ![nf-core/nascent](docs/images/nf-core-nascent_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nascent/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.157735234-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.157735234)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/nascent)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nascent-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nascent)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/nascent** is a bioinformatics best-practice analysis pipeline for nascent transcript (NT) and Transcriptional Start Site (TSS) assays.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
<!-- https://github.com/nf-core/nascent/issues/66 -->

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

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.6`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/nascent -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run nf-core/nascent --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/nascent pipeline comes with documentation about the pipeline [usage](https://nf-co.re/nascent/usage), [parameters](https://nf-co.re/nascent/parameters) and [output](https://nf-co.re/nascent/output).

## Credits

nf-core/nascent was originally written by Ignacio Tripodi ([@ignaciot](https://github.com/ignaciot)) and Margaret Gruca ([@magruca](https://github.com/magruca)).

The pipeline was re-written in Nextflow DSL2 by Edmund Miller ([@Emiller88](https://github.com/emiller88)) and Sruthi Suresh ([@sruthipsuresh](https://github.com/sruthipsuresh)) from [The Functional Genomics Laboratory](https://taehoonkim.org/) at [The Univeristy of Texas at Dallas](https://www.utdallas.edu/)

We thank the following people for their extensive assistance in the development of this pipeline:

[@apeltzer](https://github.com/apeltzer)
[@ewels](https://github.com/ewels)
[@drpatelh](https://github.com/drpatelh)
[@pditommaso](https://github.com/pditommaso)
[@FriederikeHanssen](https://github.com/FriederikeHanssen)
[Tae Hoon Kim](https://github.com/taehoonkim-phd)
[@easterwoods](https://github.com/easterwoods)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nascent` channel](https://nfcore.slack.com/channels/nascent) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/nascent for your analysis, please cite it using the following doi: [10.5281/zenodo.157735234](https://doi.org/10.5281/zenodo.157735234)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
