# TODO Plug this in
rule GM19_merge_sample_peaks:
    input:
        expand("results/2019-01-28/GM/{unit}_groseq_peak.gtf",unit=GM_SAMPLES),
    output:
        "results/2019-01-31/GM/GM19_merged_peaks.gtf"
    singularity:
        "docker://emiller88/homer:latest"
    params:
        maxDistance="100"
    shell:
        "mergePeaks -d {params.maxDistance} {input} > {output}"

# FIXME uses unmerged regions of eRNAs
rule eRNA_saf:
    input:
        "results/2018-12-02/{sample}_eRNA.bed",
    output:
        "results/2019-06-07/{sample}_eRNA.saf",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/eRNA_saf_{sample}.log"
    script:
        "../scripts/bed2saf.py"
