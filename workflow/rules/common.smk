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

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
# validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
# validate(units, schema="../schemas/units.schema.yaml")

CELLS=list(set(samples["cell"]))
GM_SAMPLES=samples[samples["cell"]=="GM"].index.tolist()
IMR_SAMPLES=samples[samples["cell"]=="IMR"].index.tolist()
TIMES=list(set(samples["time"]))
SAMPLES=list(samples["sample"])

report: "../report/workflow.rst"

# ##### wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"])

####### helpers ###########

def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]
