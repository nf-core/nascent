"""Snakemake wrapper for trimming paired-end reads using cutadapt."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


# Run command.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "featureCounts"
    " {snakemake.params.extra}"
    " -a {snakemake.params.annotation}"
    " -o {snakemake.output.counts}"
    " -T {snakemake.threads}"
    " {snakemake.input.bam} {log}")

# Move summary to expected location.
summary_path = snakemake.output.counts + '.summary'

if summary_path != snakemake.output.summary:
    shell("mv {summary_path} {snakemake.output.summary}")

