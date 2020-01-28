rule eRNA_link_genes:
    """
    List of linked eRNAs that are within 200Kb to genes
    TODO Move parameter to config
    """
    input:
        dges="results/2019-06-26/dge/rg/{cell}_de_genes.gtf",
        eRNA="results/2018-12-02/{genome}/{cell}_eRNA.bed",
    output:
        "results/2019-08-26/{genome}/{cell}_link_eRNA.bed",
    conda:
        "../../envs/bedtools.yaml"
    params:
        window=config["E-P-distance"]
    threads: 2
    shell:
        "bedtools window -u -w {params.window} -a {input.eRNA} -b {input.dges} > {output}"

rule eRNA_link_merge:
    """
    Creates a meta of eRNAs, incase transcripts differ sligtly across cell lines
    """
    input:
        GM="results/2019-08-26/hg19/GM_link_eRNA.bed",
        IMR="results/2019-08-26/hg19/IMR_link_eRNA.bed",
    output:
        "results/2019-08-26/eRNA_merged.bed",
    conda:
        "../../envs/bedtools.yaml"
    params:
        # dis="-d 2000",
        col="-c 4,5,6 -o distinct",
    shell:
        "cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin -s {params.col} > {output}"

rule eRNA_link_overlap:
    """
    Finds the overlap between the cell lines, which are virally induced
    """
    input:
        GM="results/2019-08-26/hg19/GM_link_eRNA.bed",
        IMR="results/2019-08-26/hg19/IMR_link_eRNA.bed",
    output:
        "results/2019-08-26/eRNA_overlap_viral.bed",
    conda:
        "../../envs/bedtools.yaml"
    params:
    shell:
        "bedtools intersect -wo -a {input.GM} -b {input.IMR} > {output}"

rule eRNA_link_area:
    """
    Creates a list of meta viral induced eRNAs
    """
    input:
        merge="results/2019-08-26/eRNA_merged.bed",
        ov="results/2019-08-26/eRNA_overlap_viral.bed",
    output:
        "results/2019-08-26/eRNA_viral.bed",
    conda:
        "../../envs/bedtools.yaml"
    params:
    shell:
        "bedtools intersect -wa -a {input.merge} -b {input.ov} > {output}"

rule gene_nearest_eRNA:
    input:
        dges="results/2019-06-26/dge/rg/{cell}_de_genes.gtf",
        eRNA="results/2019-06-26/dge/rg/{genome}/{cell}_de_ripgrep.bed",
    output:
        report("results/2019-11-14/{genome}/{cell}_nearest_eRNA.bed", category="Inducible Pairs")
    conda:
        "../../envs/bedtools.yaml"
    params:
    shell:
        "bedtools closest -d -a {input.dges}  -b {input.eRNA} > {output}"

rule average_e_p_distance:
    input:
        e_p="results/2019-11-14/{genome}/{cell}_nearest_eRNA.bed",
    output:
        report("results/2019-11-14/{genome}/{cell}_avg_ep.tsv", category="Inducible Pairs")
    conda:
        "../../envs/matplotlib.yaml"
    params:
    script:
        "../../scripts/nearest-avg.py"