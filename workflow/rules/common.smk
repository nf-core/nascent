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

CELLS=list(set(samples["cell"]))
GM_SAMPLES=samples[samples["cell"]=="GM"].index.tolist()
IMR_SAMPLES=samples[samples["cell"]=="IMR"].index.tolist()
UNITS=list(set(samples["time"]))
SAMPLES=list(samples["sample"])

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index("sample", drop=False)
units.index.names = ["sample_id"]
# validate(units, schema="../schemas/units.schema.yaml")


report: "../report/workflow.rst"

# ##### wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),
    # unit="|".join(units["unit"])

####### helpers ###########

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    return units.loc[wildcards.sample, "fq1" ]
