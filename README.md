# Snakemake workflow: IMR90&GM_eRNAs

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.com/Emiller88/IMR90.svg?token=4xxfcp3gAkNPaDFsDwgn&branch=master)](https://travis-ci.com/Emiller88/IMR90)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

-   Edmund Miller (@emiller88)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/IMR90/releases).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

#### Install homer

[Instructions](http://homer.ucsd.edu/homer/introduction/install.html)

Make sure to run the following to install the necessary genomes

```shell
perl /path-to-homer/configureHomer.pl -install hg18
perl /path-to-homer/configureHomer.pl -install hg19
```

In the future maybe the conda homer install will pick up on the genomes

#### Install Conda environment

With [conda](https://conda.io/en/latest/miniconda.html) installed move to the
workflow directory and execute the following.

```shell
conda env create -f environment.yml
```

#### Install R and Deseq2

From homers install guide

```md
Option 1: Use Anaconda/Bioconda to install R along with DESeq2 and EdgeR - see above. Recommended, particularly if you don't have super-user access.

Option 2: Depending on your Linux distribution, you can use a standard package manager to install samtools. Generally this option is not recommended because the version of R in the repositories is usually fairly old:
(Debian/Ubuntu): sudo apt-get install r-base r-base-dev
(Redhat/CentOS): sudo yum install r-base r-base-dev

Then run R to install Bioconductor/DESeq2/EdgeR (see below)
Option 3: Download and install R directly from the source: http://cran.cnr.berkeley.edu/
Follow the instructions to install R depending on your system.
If you picked option 2 or 3, now you'll need to run R to install DESeq2 and EdgeR:
Run R by typing "R". You may want to run this as super-user if installing for multiple users (i.e. "sudo R"). At the R prompt (should see a ">"), type the following commands:

> source("https://bioconductor.org/biocLite.R")
> biocLite()
> biocLite("DESeq2")
> biocLite("edgeR")
> q()

If you're having touble here, it might be because your version of R is too old. Consider using option 3 and get the latest stable version.
```

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.
