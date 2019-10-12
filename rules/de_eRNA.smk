rule eRNA_viral_tpm:
    input:
        "results/2019-06-03/eRNA/counts/{cell}_merged.txt"
    output:
        "results/2019-09-27/de/tpm/{cell}_tpm.txt",
    params:
    conda:
        "../envs/matplotlib.yaml"
    threads: 4
    script:
        "../scripts/tpm.py"

rule eRNA_viral_foldchange:
    input:
        "results/2019-09-27/de/tpm/{cell}_tpm.txt",
    output:
        "results/2019-09-27/de/foldchange/{cell}_foldchange.tsv",
    params:
        # cutoff=0.5
    conda:
        "../envs/matplotlib.yaml"
    threads: 4
    script:
        "../scripts/foldchange.py"
