from os import path

import numpy as np

################################################################################
# Functions                                                                    #
################################################################################

def normalize_counts(counts):
    """Normalizes expression counts using DESeq's median-of-ratios approach."""

    with np.errstate(divide="ignore"):
        size_factors = estimate_size_factors(counts)
        return counts / size_factors


def estimate_size_factors(counts):
    """Calculate size factors for DESeq's median-of-ratios normalization."""

    def _estimate_size_factors_col(counts, log_geo_means):
        log_counts = np.log(counts)
        mask = np.isfinite(log_geo_means) & (counts > 0)
        return np.exp(np.median((log_counts - log_geo_means)[mask]))

    log_geo_means = np.mean(np.log(counts), axis=1)
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col, axis=0,
        arr=counts, log_geo_means=log_geo_means)

    return size_factors


################################################################################
# Rules                                                                        #
################################################################################


def feature_counts_extra(wildcards):
    extra = config["feature_counts"]["extra"]
    if is_paired:
        extra += " -p"
    return extra


rule genes_feature_counts:
    """
    Runs featureCounts with refseq individually
    """
    input:
        bam="results/2018-10-04/{genome}/{sample}/Aligned.out.bam",
        # bai="bam/final/{sample}.bam.bai",
    output:
        counts="results/2019-06-03/{genome}/counts/per_sample/{sample}.txt",
        summary="results/2019-06-03/{genome}/qc/feature_counts/{sample}.txt"
    params:
        annotation=config["feature_counts"]["annotation"],
        extra='' #feature_counts_extra
    threads:
        config["feature_counts"]["threads"]
    log:
        "logs/{genome}/feature_counts/{sample}.txt"
    wrapper:
        "file:" + path.join(workflow.basedir, "wrappers/subread/feature_counts")

rule GM19_genes_merge_counts:
    """
    Merges individual gene counts for GM hg19
    HACK Generalize
    """
    input:
        expand("results/2019-06-03/hg19/counts/per_sample/{sample}.txt", sample=GM_SAMPLES)
    output:
        "results/2019-06-03/hg19/counts/GM_merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                        index_col=list(range(6)))
            for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)

rule IMR_genes_merge_counts:
    """
    Merges individual gene counts for IMR hg19
    HACK Generalize
    """
    input:
        expand("results/2019-06-03/hg19/counts/per_sample/{sample}.txt", sample=IMR_SAMPLES)
    output:
        "results/2019-06-03/hg19/counts/IMR_merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                        index_col=list(range(6)))
            for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)

rule genes_normalize_counts:
    """
    Take the log2 of gene counts
    """
    input:
        "results/2019-06-03/hg19/counts/{cell}_merged.txt"
    output:
        "results/2019-06-03/hg19/counts/{cell}_merged.log2.txt"
    run:
        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)))
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_counts.to_csv(output[0], sep="\t", index=True)

rule eRNA_feature_counts:
    """
    Runs featureCounts on  counts, individually
    """
    input:
        bam="results/2018-10-04/hg19/{sample}/Aligned.out.bam",
        annotation="results/2019-09-27/eRNA_viral.saf",
    output:
        counts="results/2019-06-03/hg19/counts/per_sample/{sample}_eRNA.txt",
        summary="results/2019-06-03/hg19/qc/feature_counts/{sample}_eRNA.txt"
    params:
        annotation="results/2019-09-27/eRNA_viral.saf",
        extra='-F "SAF"' #feature_counts_extra
    threads:
        config["feature_counts"]["threads"]
    log:
        "logs/{sample}/feature_counts/{sample}.txt"
    wrapper:
        "file:" + path.join(workflow.basedir, "wrappers/subread/feature_counts")

rule GM19_eRNA_merge_counts:
    """
    Merges individual eRNA counts for GM hg19
    HACK Generalize
    """
    input:
        expand("results/2019-06-03/hg19/counts/per_sample/{sample}_eRNA.txt", sample=GM_SAMPLES)
    output:
        "results/2019-06-03/hg19/counts/GM_eRNA_merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                        index_col=list(range(6)))
            for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)

rule IMR_eRNA_merge_counts:
    """
    Merges individual eRNA counts for IMR hg19
    HACK Generalize
    """
    input:
        expand("results/2019-06-03/hg19/counts/per_sample/{sample}_eRNA.txt", sample=IMR_SAMPLES)
    output:
        "results/2019-06-03/hg19/counts/IMR_eRNA_merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                        index_col=list(range(6)))
            for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)

rule eRNA_normalize_counts:
    """
    Take the log2 of eRNA counts
    """
    input:
        "results/2019-06-03/hg19/counts/{cell}_eRNA_merged.txt"
    output:
        "results/2019-06-03/hg19/counts/{cell}_eRNA_merged.log2.txt"
    run:
        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)))
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_counts.to_csv(output[0], sep="\t", index=True)
