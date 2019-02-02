# TODO Plug this in
rule GM19_merge_sample_peaks:
    input:
        expand("results/2019-01-28/GM/{unit}_groseq_peak.gtf",unit=GM_SAMPLES),
    output:
        "results/2019-01-31/GM/GM19_merged_peaks.gtf"
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    params:
        maxDistance="100"
    shell:
        "mergePeaks -d {params.maxDistance} {input} > {output}"

# FIXME uses unmerged regions of eRNAs
rule GM19_eRNA_peak:
    input:
        "results/2018-11-29/GM19_eRNA.bed",
    output:
        "results/2019-01-31/GM19_eRNA_peak.gtf",
    # FIXME Genomes don't work in a conda environment
    # conda:
    #     "../envs/homer.yaml"
    shell:
        "bed2pos.pl {input} > {output}"
