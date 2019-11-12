rule eRNA_gene_groups:
    """
    Report groups of Genes and their possible eRNAs
    FIXME Cluttered by multiple exons
    """
    input:
        dges="results/2019-06-26/dge/rg/{cell}_de_genes.gtf",
        eRNA="results/2019-06-26/dge/rg/{genome}/{cell}_de_ripgrep.bed",
    output:
        report("results/2019-08-26/{genome}/{cell}_eRNA_gene_group.bed", category="Inducible Pairs")
    conda:
        "../../envs/bedtools.yaml"
    params:
        window="200000"
    threads: 2
    shell:
        "bedtools window -w {params.window} -a {input.dges}  -b {input.eRNA} > {output}"

rule group_names:
    """
    Creates a header
    Prints the gene_id and eRNA_id
    Cleans the gene_id
    Changes space to tab
    FIXME Removes none unique genes and eRNA groups
    """
    input:
        "results/2019-08-26/{genome}/{cell}_eRNA_gene_group.bed",
    output:
        report("results/2019-10-29/{genome}/{cell}_group_names.tsv", category="Inducible Pairs")
    conda:
        "../../envs/gawk.yaml"
    shell:
        """echo -e "Gene_id \teRNA_id" > {output} &&\
        awk -F '\\t' '{{print $9, $13}}' {input} \
        | sed 's/[";]//g' - \
        | awk -F " " '{{OFS="\\t";}} {{print $2, $NF}}' - \
        | uniq - >> {output}"""

rule common_gene_group:
    input:
        groups=expand("results/2019-10-29/{genome}/{cell}_group_names.tsv", genome=["hg19"],cell=CELLS),
    output:
        report("results/2019-10-29/{genome}/common_group_names.tsv", category="Inducible Pairs")
    conda:
        "../../envs/gawk.yaml"
    shell:
        """awk 'NR==FNR{{a[$1]=1;next}}a[$1]' {input.groups} > {output}"""
