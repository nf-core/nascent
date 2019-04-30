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
rule GM19_eRNA_gff:
    input:
        "results/2018-11-29/GM19_eRNA.bed",
    output:
        "results/2019-01-31/GM19_eRNA_peak.gff",
    conda:
        "../../envs/gawk.yaml"
    shell:
        # https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html
        """awk '{{print $1"\\thomer\\tenhancer\\t"($2+1)"\\t"$3"\
        \\t"$5"\\t"$6"\\t"$9"\\t"(substr($0, index($0,$10)))\
        }}' {input} > {output}"""
