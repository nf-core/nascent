"""
eRNA_GM_hg19(legacy eRNAs) is identified with hg18(using hg18 uniqmap) -> liftOver -> hg19
"""
rule test_eRNA_vs_Peng:
    """
    Compares identified eRNAs to legacy eRNAs
    """
    input:
        edmund="results/2018-12-02/{genome}/{cell}_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
    output:
        report("results/2018-11-10/test/{genome}/{cell}_eRNA_overlaps.bed", category="eRNA Prediction")
    log:
        "logs/{genome}/test_{cell}_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.edmund} \
        -sorted -u > {output} 2> {log}"

rule test_GM18_liftOver_vs_Peng:
    """
    Compares GM18 eRNAs that are lifted over to hg19 to legacy eRNAs
    """
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
    output:
        report("results/2018-11-10/test/hg19_eRNA_overlaps.bed", category="eRNA Prediction")
    log:
        "logs/hg18/test_GM18_liftOver_eRNAs.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.peng} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_hg19_vs_hg18:
    """
    Compares indentified hg18 to hg19 in one cell line
    """
    input:
        hg19="results/2018-12-02/hg19/{cell}_eRNA.bed",
        hg18="results/2018-12-02/hg18/{cell}_eRNA.bed",
    output:
        report("results/2018-11-10/test/{cell}_hg19_vs_hg18_eRNA.bed", category="eRNA Prediction")
    log:
        "logs/hg19/test_{cell}_eRNAs_hg19_vs_hg18.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.hg19} -b {input.hg18} \
        -sorted -u > {output} 2> {log}"

rule test_eRNA_vs_liftOver:
    """
    Test identified eRNAs to liftedOver GM
    TODO Generalize this to create liftOver for both cell lines
    """
    input:
        liftOver="results/2018-11-10/eRNA_GM_liftover_hg19.sorted.bed",
        eRNA="results/2018-12-02/{genome}/{cell}_eRNA.bed",
    output:
        report("results/2018-11-10/test/{genome}/{cell}_eRNA_vs_liftOver.bed", category="eRNA Prediction")
    log:
        "logs/{genome}/test_{cell}_eRNAs_liftover.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.eRNA} -b {input.liftOver} \
        -sorted -u > {output} 2> {log}"

rule test_IMR_vs_GM:
    """
    Compares identified eRNAs across cell lines
    """
    input:
        IMR="results/2018-12-02/{genome}/IMR_eRNA.bed",
        GM19="results/2018-12-02/{genome}/GM_eRNA.bed",
    output:
        report("results/2018-11-10/test/{genome}/IMR_eRNA_vs_GM_{genome}.bed", category="eRNA Prediction")
    log:
        "logs/{genome}/test_IMR_vs_GM.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.GM19} \
        -sorted -u > {output} 2> {log}"

rule fig_predicted_eRNA_peng:
    """
    Creates a venn diagram of eRNA transcripts from rule test_eRNA_vs_Peng
    """
    input:
        edmund="results/2018-12-02/{genome}/{cell}_eRNA.bed",
        peng="data/2018-01-25/eRNA_GM_hg19.bed",
        overlap="results/2018-11-10/test/{genome}/{cell}_eRNA_overlaps.bed",
    output:
        report("results/2018-10-12/{genome}/{cell}_eRNA_overlaps.svg", category="Figures")
    log:
        "logs/{genome}/figure/{cell}_eRNAs_vs_peng.log"
    conda:
        "../../envs/venn.yaml"
    params:
        title="Predicted eRNAs vs Peng's GM19",
    script:
        "../../scripts/venn-smk.py"

rule fig_predicted_eRNA_cross_cell:
    """
    Creates a venn diagram of eRNA transcripts across cell lines from rule test_IMR_vs_GM
    """
    input:
        IMR="results/2018-12-02/{genome}/IMR_eRNA.bed",
        GM19="results/2018-12-02/{genome}/GM_eRNA.bed",
        overlap="results/2018-11-10/test/{genome}/IMR_eRNA_vs_GM_{genome}.bed",
    output:
        report("results/2018-10-12/{genome}/eRNA_cross_cell.svg", category="Figures")
    log:
        "logs/{genome}/figure/eRNA_cross_cell.log"
    conda:
        "../../envs/venn.yaml"
    params:
        title="IMR vs GM",
    script:
        "../../scripts/venn-smk.py"

rule test_IMR_vs_GM_0h:
    """
    Compares identified eRNAs across cell lines
    """
    input:
        IMR="results/2018-12-02/{genome}/sample_IMR0h_IMR_eRNA.bed",
        GM="results/2018-12-02/{genome}/sample_GM0h_GM_eRNA.bed",
    output:
        report("results/2018-11-10/test/{genome}/IMR_eRNA_vs_GM_0h_{genome}.bed", category="eRNA Prediction")
    log:
        "logs/{genome}/test_IMR_vs_GM_0h.log"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.IMR} -b {input.GM} \
        -sorted -u > {output} 2> {log}"

rule fig_predicted_eRNA_cross_cell_0h:
    """
    Creates a venn diagram of eRNA transcripts across cell lines from rule test_IMR_vs_GM
    """
    input:
        IMR="results/2018-12-02/{genome}/sample_IMR0h_IMR_eRNA.bed",
        GM="results/2018-12-02/{genome}/sample_GM0h_GM_eRNA.bed",
        overlap="results/2018-11-10/test/{genome}/IMR_eRNA_vs_GM_0h_{genome}.bed",
    output:
        report("results/2018-10-12/{genome}/eRNA_cross_cell_0h.svg", category="Figures")
    log:
        "logs/{genome}/figure/eRNA_cross_cell.log"
    conda:
        "../../envs/venn.yaml"
    params:
        title="IMR 0h vs GM 0h",
    script:
        "../../scripts/venn-smk.py"
