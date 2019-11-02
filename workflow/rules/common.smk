import pandas as pd
from snakemake.utils import validate, min_version
from itertools import product
##### set minimum snakemake version #####
min_version("5.5.0")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3:4.6.14"

##### load config and sample sheets #####
configfile: "config/config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype=str, sep='\t').set_index(["cell", "name"], drop=False,)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
validate(samples, schema="schemas/samples.schema.yaml")

CELLS=list(set(samples["cell"]))
GM_SAMPLES=samples.loc[("GM")].index.tolist()
IMR_SAMPLES=samples.loc[("IMR")].index.tolist()
UNITS=list(set(samples["time"]))
SAMPLES=list(samples["name"])

report: "report/workflow.rst"

##### wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"])

####### helpers ###########

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]
