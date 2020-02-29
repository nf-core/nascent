# ![nf-core/nascent](docs/images/nf-core-nascent_logo.png)

**Nascent Transcription Processing Pipeline**.

[![GitHub Actions CI Status](https://github.com/nf-core/nascent/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/nascent/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/nascent/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/nascent/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nascent.svg)](https://hub.docker.com/r/nfcore/nascent)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/nascent -profile test,<docker/singularity/conda/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

<!-- TODO nf-core: Update the default command above used to run the pipeline -->

```bash
nextflow run nf-core/nascent -profile <docker/singularity/conda/institute> --reads '*_R{1,2}.fastq.gz' --genome GRCh37
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The nf-core/nascent pipeline comes with documentation about the pipeline, found in the `docs/` directory:

If you have used this pipeline in your research, please cite it using the DOI mentioned above.

The nf-core/nascent pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

This pipeline is designed to process the sequencing output of nascent transcription assays, like GRO-seq or PRO-seq. It produces bedGraph- and bigWig-fomatted outputs after mapping strand-specific reads, as well as other useful outputs like quality control reports or IGV-ready (Integrative Genomics Viewer) tdf files.

### Quick start

Edit the appropriate config file, e.g. `conf/slurm_grch38.config`, to ensure the proper paths are set for genome reference files and other executables (look for all mentions of `COMPLETE_*`). Variable names should hopefully be self-explanatory. You can specify the Nextflow working directory and output directory with flags. Note you must also now specify the email to which the report will be sent for the run.

```console
nextflow run nf-core/nascent --reads '*_R{1,2}.fastq.gz' -profile standard,docker
```

## Arguments

### Required Arguments

| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                       |
| --fastqs  | \</project/\*\_{R1,R2}\*.fastq\> | Directory pattern for fastq files.                                   |
| --sras    | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --genome_id | \<'hg38'>                      | Genome ID to which the samples will be mapped (e.g. hg38, mm10, rn6).|
| --workdir | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
| --email   | \<EMAIL\>                        | Where to send workflow report email.                                 |

### Save Options

| Arguments  | Usage         | Description                                               |
|------------|---------------|-----------------------------------------------------------|
| --outdir   | \</project/\> | Specifies where to save the output from the nextflow run. |
| --savefq   |               | Compresses and saves raw fastq reads.                     |
| --saveTrim |               | Compresses and saves trimmed fastq reads.                 |
| --saveAll  |               | Compresses and saves all fastq reads.                     |
| --skipBAM  |               | Skips saving BAM files (only save CRAM). Default=False    |

### Input File Options

| Arguments    | Usage       | Description                                                                  |
|--------------|-------------|------------------------------------------------------------------------------|
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |
| --flip       |             | Reverse complements each strand. Necessary for some library preps.           |

### Performance Options

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --threadfqdump  |             | Runs multi-threading for fastq-dump for sra processing. |

### QC Options

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC.                                   |
| --skipRSeQC     |             | Skip running RSeQC.                                     |

## Credits

nf-core/nascent was originally written by Ignacio Tripodi ([@ignaciot](https://github.com/ignaciot)) and Margaret Gruca ([@magruca](https://github.com/magruca)).

Many thanks to the nf-core team and all who provided invaluable feedback and assistance along the way, particularly to [@apeltzer](https://github.com/apeltzer), [@ewels](https://github.com/ewels), [@drpatelh](https://github.com/drpatelh), and [@pditommaso](https://github.com/pditommaso).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/nascent) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/nascent for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
