import pandas as pd
from snakemake.utils import validate
from itertools import product

# singularity: "docker://continuumio/miniconda3:4.6.14"

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], dtype=str).set_index(["cell", "name"], drop=False)
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
        # fastqc
        fastqc_zip=expand("results/2018-10-03/fastqc/{cell}{unit}.zip", filtered_product, cell=CELLS, unit=UNITS,),
        # GM hg18
        # eRNAliftOverhg19="results/2018-11-25/hg19_eRNA_overlaps.bed",
        # GM18_test1="results/2018-11-10/hg18_eRNA_overlaps.bed",
        # GM18_test2="results/2018-11-25/hg19_eRNA_overlaps.bed",
        # GM hg19
        GM_counts=expand("results/2019-01-28/GM/{unit}_groseq_peak.bed", unit=GM_SAMPLES),
        # hg19eRNA="results/2018-11-29/eRNA_GM_hg19.bed",
        # GM19_test1="results/2018-11-30/hg19_eRNA_overlaps_Peng.bed",
        # GM19_test2="results/2018-11-30/hg19_eRNA_overlaps_liftOver.bed",
        # IMR hg19
        IMR_counts=expand("results/2019-01-28/IMR/{unit}_groseq_peak.bed", unit=IMR_SAMPLES),
        # IMR_eRNA="results/2018-12-02/eRNA_IMR_hg19.bed",
        IMR_test1="results/2018-11-10/test/IMR_eRNA_vs_Peng.bed",
        IMR_test2="results/2018-11-10/test/IMR_eRNA_vs_GM19.bed",
        IMR_test3="results/2018-11-10/test/IMR_eRNA_vs_liftOver.bed",
        # Differential Analysis
        # GM_diff="results/2018-01-30/GM19_eRNA_diffPeaks.txt",
        # IMR_diff="results/2018-01-30/IMR_eRNA_diffPeaks.txt",
        # GM_annotation=expand("results/2019-02-05/GM/{unit}_outputannotation.txt", unit=GM_SAMPLES),
        # genes_limma=expand("results/2019-06-26/dge/limma/{cell}_{fig}_limma.png", cell=["GM19", "IMR"], fig=["fig1", "fig2", "fig3", ]),
        # erna_limma=expand("results/2019-06-26/eRNA/limma/{cell}_{fig}_limma.png", cell=["GM19", "IMR"], fig=["fig1", "fig2", "fig3", ]),
        genes_foldchange=expand("results/2019-06-26/dge/foldchange/{cell}_foldchange.tsv", cell=["GM19", "IMR"]),
        linkedeRNAs=expand("results/2019-08-26/{cell}_link_eRNA.bed", cell = ["GM19", "IMR"]),
        merge="results/2019-08-26/eRNA_viral.bed",
        eRNAcounts=expand("results/2019-06-03/eRNA/counts/{cell}_merged.txt", cell = ["GM19", "IMR"]),
        eRNA_foldchange=expand("results/2019-09-27/de/foldchange/{cell}_foldchange.tsv", cell=["GM19", "IMR"]),
        l2=expand("results/2019-08-26/{cell}_L2.bed", cell=["GM19","IMR"]),
        genesl2=expand("results/2019-10-01/{cell}_genes_L2.bed", cell=["GM19","IMR"]),
        baml2=expand("results/2019-10-06/l2/{sample}_bam_L2.bed", sample=IMR_SAMPLES),
        l2report="results/2019-10-01/L2.tsv",

include: "rules/data.smk"

include: "rules/GM18/bowtie2.smk"
include: "rules/GM18/homer.smk"
include: "rules/GM18/removeGenes.smk"
include: "rules/GM18/keepHistones.smk"
include: "rules/GM18/liftOver.smk"
include: "rules/GM18/test_peng_eRNAs.smk"

include: "rules/GM19/fastqc.smk"
include: "rules/GM19/bowtie2.smk"
include: "rules/GM19/homer.smk"
include: "rules/GM19/removeGenes.smk"
include: "rules/GM19/keepHistones.smk"
include: "rules/GM19/test_eRNAs.smk"

include: "rules/IMR/fastqc.smk"
include: "rules/IMR/bowtie2.smk"
include: "rules/IMR/homer.smk"
include: "rules/IMR/removeGenes.smk"
include: "rules/IMR/keepHistones.smk"
include: "rules/IMR/test_eRNAs.smk"

include: "rules/fastqc.smk"

include: "rules/eRNAcleaning.smk"
include: "rules/countReads.smk"
include: "rules/dge.smk"
include: "rules/linkedRNAs.smk"
include: "rules/de_eRNA.smk"
include: "rules/test_eRNA.smk"

include: "rules/reports/eRNA.smk"
