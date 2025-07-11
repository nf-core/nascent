# nf-core/nascent: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/nascent/usage](https://nf-co.re/nascent/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

:::info
The sample column is essentially a concatenation of the group and replicate columns. If all values of sample have the same number of underscores, fields defined by these underscore-separated names may be used in the transcript identification produced by the pipeline, to regain the ability to represent different groupings.

`GM_0h` and `GM_1h` would be grouped for example but `GM0h` and `GM1h` would go through individual transcript identification
:::

## Alignment Options

By default, the pipeline uses [BWA](https://bio-bwa.sourceforge.net/) (i.e. `--aligner bwa`) to map the raw FastQ reads to the reference genome. Research as to which aligner works best with Nascent Transcript and Transcription Start Site assays is pending.

## Reference genome files

The minimum reference genome requirements are a FASTA and GTF file, all other files required to run the pipeline can be generated from these files. However, it is more storage and compute friendly if you are able to re-use reference genome files as efficiently as possible. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices (e.g. those unavailable on [AWS iGenomes](https://nf-co.re/usage/reference_genomes)) so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space. You can then either provide the appropriate reference genome files on the command-line via the appropriate parameters (e.g. `--star_index '/path/to/BWA/index/'`) or via a custom config file.

- If `--genome` is provided then the FASTA and GTF files (and existing indices) will be automatically obtained from AWS-iGenomes unless these have already been downloaded locally in the path specified by `--igenomes_base`.
- If `--gff` is provided as input then this will be converted to a GTF file, or the latter will be used if both are provided.

## Quantification Options

Currently only featureCounts is supported for quantification. It counts both the genes, and the predicted transcripts.

## Transcript Identification Options

The current options for transcript identification include [GroHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html), [HOMER](http://homer.ucsd.edu/), and [PINTS](https://pints.yulab.org/).

The default transcript identification option is PINTS, and HOMER if the transcript `assay_type` is `GROseq` but this may change in future releases.

### PINTS

PINTS handles the majority of the transcript identification, since it covers all of the supported assays.

PINTS can use a lot of memory while running, so [a scatter-gather pattern was implemented](https://github.com/nf-core/nascent/blob/136a9ca2390121639e823e39e508afe9b6970d77/subworkflows/local/transcript_identification/main.nf#L47-L74).

It splits the identification up by the chromosomes available in the provided FASTA file. Some of the chromosomes are skipped because PINTS throws an error when it doesn't find any regions. If this causes an issue with your analysis please open an issue.

### GroHMM

groHMM is split into two steps: parameter tuning and transcript identification.

When running the pipeline with groHMM as a transcript identification method, the pipeline will automatically perform a parameter tuning process. This process is unique to the groHMM transcript identification method and is designed to select the optimal hold-out parameters for the groHMM algorithm. See [this issue](https://github.com/dankoc/groHMM/issues/4) for more information.

In the groHMM vignette, the code is ran using a single mclapply call, which is a scatter gather approach. This is not ideal for large datasets, because it ends up being bottle-necked by the memory available on your local machine. To improve this, we have written a Nextflow script that runs the pipeline with a scatter gather approach. This is done by running the pipeline with a single hold-out parameter, and then the next parameter, and so on. This is more memory efficient and scales better to larger datasets. The results are then combined in the end as intended and used in the transcript identification process.

#### groHMM Parameters

> The detectTranscripts function also uses two hold-out parameters. These parameters, specified by the arguments LtProbB and UTS, represents the log-transformed transition probability of switching from transcribed state to non-transcribed state and variance of the emission probability for reads in the non-transcribed state, respectively. Holdout parameters are used to optimize the performance of HMM predictions on known genes.

In the pipeline, the parameters are specified as follows:
grohmm_min_uts = 5
grohmm_max_uts = 45
grohmm_min_ltprobb = -100
grohmm_max_ltprobb = -400

Which will then create a job for each parameter combination. For example (5,-100), (5,-150), (10,-100), (10,-150)...

If you have indentified a good set of parameters, you can run the pipeline with those parameters by specifying, all 4 values.

For example if you have indentified that the best parameters for your data are 15,-200:

```json
{
  "grohmm_min_uts": 15,
  "grohmm_max_uts": 15,
  "grohmm_min_ltprobb": -200,
  "grohmm_max_ltprobb": -200
}
```

### Homer

HOMER is used for transcript identification when the `assay_type` is set to `GROseq`. HOMER's GRO-seq analysis capabilities include:

- De novo transcript identification from GRO-seq data
- Support for analyzing nascent RNA production
- Detection of various RNA species including:
  - Protein coding transcripts
  - Promoter anti-sense transcripts
  - Enhancer templated transcripts (eRNAs)
  - Long and short non-coding RNAs
  - miRNA transcripts
  - Pol III and Pol I transcripts

HOMER uses uniquely mappable regions to improve transcript detection in repetitive regions. The pipeline can automatically download the appropriate uniqmap files for supported genomes:

- Human: hg19, hg38
- Mouse: mm10
- Fly: dm6

:::info
**This setting is off by default**
:::

To find the full list of uniqmaps supplied by the author check http://homer.ucsd.edu/homer/data/uniqmap/. To build a uniqmap for a genome that isn't supported, check out [homer-uniqmap-nf](https://github.com/Functional-Genomics-Lab/homer-uniqmap-nf).

The transcript detection algorithm:

1. Tracks along each strand looking for continuous GRO-seq signal
2. Starts transcripts when encountering high read density
3. Stops transcripts when signal decreases significantly
4. Creates new transcripts when signal increases sustainably
5. Filters out artifactual spikes that don't extend over distance

Key parameters that can be tuned:

- tssFold: Fold change required at transcript start (default: 4)
- bodyFold: Fold change required in transcript body (default: 3)
- minBodySize: Minimum transcript body size (default: 600bp)
- maxBodySize: Maximum transcript body size (default: 10000bp)

```nextflow
withName: HOMER_FINDPEAKS {
    ext.args = "-style groseq -tssFold 4 -bodyFold 3"
}
```

For more info check the [Homer GRO-seq Tutorial](http://homer.ucsd.edu/homer/ngs/groseq/groseq.html).

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/nascent --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/nascent -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/nascent
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/nascent releases page](https://github.com/nf-core/nascent/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
