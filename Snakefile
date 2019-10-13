import pandas as pd
from snakemake.utils import validate, min_version
from itertools import product
##### set minimum snakemake version #####
min_version("5.5.0")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype=str, sep='\t').set_index(["cell", "name"], drop=False,)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
validate(samples, schema="schemas/samples.schema.yaml")

CELLS=list(set(samples["cell"]))
GM_SAMPLES=samples.loc[("GM")].index.tolist()
IMR_SAMPLES=samples.loc[("IMR")].index.tolist()
UNITS=list(set(samples["time"]))
SAMPLES=list(samples["name"])

def filter_combinator(combinator, blacklist):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            if frozenset(wc_comb) not in blacklist:
                yield wc_comb
    return filtered_combinator

forbidden = {frozenset(wc_comb.items()) for wc_comb in config["missing"]}
filtered_product = filter_combinator(product, forbidden)

report: "report/workflow.rst"

rule all:
    input:
        aws="scripts/aws-igenomes.sh",
        # fastqc
        fastqc_zip=expand("results/2018-10-03/fastqc/{cell}{unit}.zip", filtered_product, cell=CELLS, unit=UNITS,),
        # eRNA Prediction
        eRNApeng=expand("results/2018-11-10/test/{genome}/{cell}_eRNA_overlaps.bed", genome=config["genomes"], cell=CELLS,),
        liftoverPeng="results/2018-11-10/test/hg19_eRNA_overlaps.bed",
        hg19vshg18=expand("results/2018-11-10/test/{cell}_hg19_vs_hg18_eRNA.bed", cell=CELLS,),
        predictionVsLiftOver=expand("results/2018-11-10/test/{genome}/{cell}_eRNA_vs_liftOver.bed", genome=config["genomes"], cell=CELLS,),
        IMRvsGM="results/2018-11-10/test/IMR_eRNA_vs_GM19.bed",
        # Figures
        VennvsPeng=expand("results/2018-10-12/{genome}/{cell}_eRNA_overlaps.svg", genome=config["genomes"], cell=CELLS,),
        # Differential Analysis
        # GM_diff="results/2018-01-30/GM19_eRNA_diffPeaks.txt",
        # IMR_diff="results/2018-01-30/IMR_eRNA_diffPeaks.txt",
        # GM_annotation=expand("results/2019-02-05/GM/{unit}_outputannotation.txt", unit=GM_SAMPLES),
        # genes_limma=expand("results/2019-06-26/dge/limma/{cell}_{fig}_limma.png", cell=["GM19", "IMR"], fig=["fig1", "fig2", "fig3", ]),
        # erna_limma=expand("results/2019-06-26/eRNA/limma/{cell}_{fig}_limma.png", cell=["GM19", "IMR"], fig=["fig1", "fig2", "fig3", ]),
        genes_foldchange=expand("results/2019-06-26/dge/foldchange/{cell}_foldchange.tsv", cell=["GM", "IMR"]),
        linkedeRNAs=expand("results/2019-08-26/hg19/{cell}_link_eRNA.bed", cell = ["GM", "IMR"]),
        merge="results/2019-08-26/eRNA_viral.bed",
        eRNAcounts=expand("results/2019-06-03/eRNA/counts/{cell}_merged.txt", cell = ["GM", "IMR"]),
        eRNA_foldchange=expand("results/2019-09-27/de/foldchange/{cell}_foldchange.tsv", cell=["GM", "IMR"]),
        l2=expand("results/2019-08-26/{cell}_L2.bed", cell=["GM","IMR"]),
        genesl2=expand("results/2019-10-01/{cell}_genes_L2.bed", cell=["GM","IMR"]),
        l2report="results/2019-10-01/L2.tsv",

include: "rules/data.smk"

###############
# Reproducing #
###############
include: "rules/GM18/bowtie2.smk"
include: "rules/GM18/homer.smk"
include: "rules/GM18/liftOver.smk"

###################
# eRNA Prediction #
###################
include: "rules/fastqc.smk"
include: "rules/bowtie2.smk"
include: "rules/homer.smk"
include: "rules/removeGenes.smk"
include: "rules/keepHistones.smk"
include: "rules/test_eRNA_prediction.smk"

###################
# Inducible Pairs #
###################
include: "rules/eRNAcleaning.smk"
include: "rules/countReads.smk"
include: "rules/dge.smk"
include: "rules/linkedRNAs.smk"
include: "rules/de_eRNA.smk"
include: "rules/test_eRNA.smk"
