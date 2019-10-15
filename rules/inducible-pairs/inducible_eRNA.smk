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

rule eRNA_Inducible_id:
    input:
        inducible="results/2019-09-27/de/foldchange/{cell}_foldchange.tsv",
    output:
        "results/2019-08-26/dge/rg/{cell}_de_eRNA_id.txt",
    conda:
        "../envs/gawk.yaml"
    threads: 1
    shell:
        """awk -F \"\\t\" '{{ if (NR!=1){{ print $1 }}}}' {input} > {output}"""

rule eRNA_ripgrep_id:
    input:
        inducibleId="results/2019-08-26/dge/rg/{cell}_de_eRNA_id.txt",
        linkedeRNAs="results/2019-08-26/{genome}/{cell}_link_eRNA.bed",
    output:
        "results/2019-06-26/dge/rg/{genome}/{cell}_de_ripgrep.bed",
    conda:
        "../envs/ripgrep.yaml"
    threads: 4
    shell:
        """rg --dfa-size-limit 2G -w -f {input.geneid} {input.refseq} > {output}"""
