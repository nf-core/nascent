# Snakemake workflow: IMR90&GM_eRNAs

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
![Tests](https://github.com/Emiller88/eRNA-GRO-Seq/workflows/Tests/badge.svg)
[![pipeline status](https://gitlab.com/functional-genomics/analysis/eRNA-GRO-Seq/badges/develop/pipeline.svg)](https://gitlab.com/functional-genomics/analysis/eRNA-GRO-Seq/commits/develop)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

-   Edmund Miller (@emiller88)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/emiller88/eRNA-GRO-Seq/releases).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

#### Install homer(Optional)

Using `--use-singularity` will download the docker container with the genomes
installed, it's a huge containter however.

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

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config/config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake -j 999 --use-conda --cluster-config config/cluster.json \
        --cluster "sbatch -A {cluster.account} -p {cluster.partition}\
        -n {cluster.n}  -t {cluster.time}"

or

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.
