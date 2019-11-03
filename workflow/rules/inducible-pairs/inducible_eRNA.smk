rule eRNA_viral_tpm:
    """
    Creates tpm values
    HACK Generalize with genes
    """
    input:
        "results/2019-06-03/hg19/counts/{cell}_eRNA_merged.txt"
    output:
        "results/2019-09-27/de/tpm/{cell}_eRNA_tpm.txt",
    params:
    conda:
        "../../envs/matplotlib.yaml"
    threads: 4
    script:
        "../../scripts/tpm.py"

rule inducible_eRNA_viral:
    """
    Determines inducible genes with a -1 < foldchange > 1 & under < 24h
    TODO Add one-step autocorrelation
    TODO Bring in settings from config
    """
    input:
        "results/2019-09-27/de/tpm/{cell}_eRNA_tpm.txt",
    output:
        "results/2019-09-27/de/foldchange/{cell}_eRNA_foldchange.tsv",
    params:
        cutoff=config["inducible"]["foldchange"]
    conda:
        "../../envs/matplotlib.yaml"
    threads: 4
    script:
        "../../scripts/foldchange.py"

rule eRNA_Inducible_id:
    """
    Creates a list of eRNA ids
    """
    input:
        inducible="results/2019-09-27/de/foldchange/{cell}_eRNA_foldchange.tsv",
    output:
        "results/2019-08-26/dge/rg/{cell}_de_eRNA_id.txt",
    conda:
        "../../envs/gawk.yaml"
    threads: 1
    shell:
        """awk -F \"\\t\" '{{ if (NR!=1){{ print $1 }}}}' {input} > {output}"""

rule eRNA_ripgrep_id:
    """
    Finds eRNAs BED entries that are inducible
    """
    input:
        inducibleId="results/2019-08-26/dge/rg/{cell}_de_eRNA_id.txt",
        linkedeRNAs="results/2019-08-26/{genome}/{cell}_link_eRNA.bed",
    output:
        "results/2019-06-26/dge/rg/{genome}/{cell}_de_ripgrep.bed",
    conda:
        "../../envs/ripgrep.yaml"
    threads: 4
    shell:
        """rg --dfa-size-limit 2G -w -f {input.inducibleId} {input.linkedeRNAs} > {output}"""
