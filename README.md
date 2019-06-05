# ![nfcore/nascent](docs/images/nascent_logo.png)

[![Build Status](https://travis-ci.com/nf-core/nascent.svg?branch=master)](https://travis-ci.com/nf-core/nascent)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/157735234.svg)](https://zenodo.org/badge/latestdoi/157735234)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nascent.svg)](https://hub.docker.com/r/nfcore/nascent)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
This nascent transcription processing pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


#### Reference

If you've used this pipeline in your research, you can cite this pipeline using DOI 10.17605/OSF.IO/SV4UB ([OSF project](https://osf.io/sv4ub/)).

### Documentation
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

    nextflow run nf-core/nascent --reads '*_R{1,2}.fastq.gz' -profile standard,docker

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
