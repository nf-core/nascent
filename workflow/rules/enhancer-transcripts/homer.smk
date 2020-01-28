rule IMR_meta_makeTagDirectory:
    """
    Creates a tag directory using homer and checkes GC content
    HACK Generalize
    """
    input:
        expand("results/2018-10-04/hg19/{unit}/Aligned.out.sam",unit=IMR_SAMPLES)
    output:
        directory("results/2018-11-07/hg19/IMR_meta_tagDir"),
    singularity: config["homer"]["makeTagDir"]["singularity"]
    params:
        genome = "hg19",
        extra = config["homer"]["makeTagDir"]["extra"],
    shell:
        "makeTagDirectory {output} -genome {params.genome} {params.extra} {input}"

rule GM_meta_makeTagDirectory:
    """
    Creates a tag directory using homer and checkes GC content
    HACK Generalize
    """
    input:
        expand("results/2018-10-04/hg19/{unit}/Aligned.out.sam",unit=GM_SAMPLES)
    output:
        directory("results/2018-11-07/hg19/GM_meta_tagDir")
    singularity: config["homer"]["makeTagDir"]["singularity"]
    params:
        genome = "hg19",
        extra = config["homer"]["makeTagDir"]["extra"],
    shell:
        "makeTagDirectory {output} -genome {params.genome} {params.extra} {input}"

rule meta_makeUCSCfile:
    """
    Creates a UCSC genome browser file using homer from a tag directory
    FIXME Broken
    """
    input:
        "results/2018-11-07/{genome}/{cell}_meta_tagDir"
    output:
        report("results/2018-11-07/{genome}/{cell}_{genome}_ucsc.gz", category="homer")
    singularity: config["homer"]["makeTagDir"]["singularity"]
    params:
        fixedOutput="results/2018-11-07/{genome}/{cell}_{genome}_ucsc",
        skips=config["homer"]["UCSC"]["skipChr"],
    shell:
        "makeUCSCfile {input} -o {params.fixedOutput} -strand separate -skipChr {params.skips}"

rule hg19_meta_findPeaks:
    """
    Uses homer findPeaks to indentify GRO-Seq transcripts
    HACK Generalize
    """
    input:
        tagdir="results/2018-11-07/hg19/{cell}_meta_tagDir",
        uniqmap=config["homer"]["findPeaks"]["uniqmap"],
    output:
        "results/2018-11-07/hg19/{cell}_meta_transcripts.txt"
    singularity: config["homer"]["findPeaks"]["singularity"]
    params:
        style = config["homer"]["findPeaks"]["style"],
        tssSize = config["homer"]["findPeaks"]["tssSize"],
        minBodySize = config["homer"]["findPeaks"]["minBodySize"],
        tssFold = config["homer"]["findPeaks"]["tssFold"],
        bodyFold = config["homer"]["findPeaks"]["bodyFold"],
        endFold = config["homer"]["findPeaks"]["endFold"],
        fragLength = config["homer"]["findPeaks"]["fragLength"],
        confPvalue = config["homer"]["findPeaks"]["confPvalue"],
        pseudoCount = config["homer"]["findPeaks"]["pseudoCount"],
    log: "logs/hg19/homer/{cell}_findPeaks"
    shell:
        "findPeaks {input.tagdir} -style {params.style} -uniqmap {input.uniqmap}\
        -tssSize {params.tssSize} -minBodySize {params.minBodySize} -tssFold {params.tssFold}\
        -bodyFold {params.bodyFold} -endFold {params.endFold} -confPvalue {params.confPvalue} -pseudoCount {params.pseudoCount}\
        -o {output}"

rule homer_meta_pos2bed:
    """
    Coverts transcripts from homer format to bed
    """
    input:
        "results/2018-11-07/{genome}/{cell}_meta_transcripts.txt"
    output:
        "results/2018-11-07/{genome}/{cell}_meta_transcripts.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} | sort -k1,1 -k2,2n - > {output}"

rule sample_makeTagDirectory:
    """
    Creates a tag directory using homer for 0h
    """
    input:
        "results/2018-10-04/{genome}/{sample}/Aligned.out.sam"
    output:
        directory("results/2018-11-07/{genome}/{sample}_tagDir"),
    singularity: config["homer"]["makeTagDir"]["singularity"]
    params:
        genome = "hg19",
        extra = config["homer"]["makeTagDir"]["extra"],
    shell:
        "makeTagDirectory {output} -genome {params.genome} {params.extra} {input}"

rule sample_findPeaks:
    """
    Uses homer findPeaks to indentify GRO-Seq transcripts
    HACK Generalize
    """
    input:
        tagdir="results/2018-11-07/{genome}/{sample}_tagDir",
        uniqmap=config["homer"]["findPeaks"]["uniqmap"],
    output:
        "results/2018-11-07/{genome}/{sample}_transcripts.txt"
    singularity: config["homer"]["findPeaks"]["singularity"]
    params:
        style = config["homer"]["findPeaks"]["style"],
        tssSize = config["homer"]["findPeaks"]["tssSize"],
        minBodySize = config["homer"]["findPeaks"]["minBodySize"],
        tssFold = config["homer"]["findPeaks"]["tssFold"],
        bodyFold = config["homer"]["findPeaks"]["bodyFold"],
        endFold = config["homer"]["findPeaks"]["endFold"],
        fragLength = config["homer"]["findPeaks"]["fragLength"],
        confPvalue = config["homer"]["findPeaks"]["confPvalue"],
        pseudoCount = config["homer"]["findPeaks"]["pseudoCount"],
    shell:
        "findPeaks {input.tagdir} -style {params.style} -uniqmap {input.uniqmap}\
        -tssSize {params.tssSize} -minBodySize {params.minBodySize} -tssFold {params.tssFold}\
        -bodyFold {params.bodyFold} -endFold {params.endFold} -confPvalue {params.confPvalue} -pseudoCount {params.pseudoCount}\
        -o {output}"

rule sample_pos2bed:
    """
    Coverts transcripts from homer format to bed
    """
    input:
        "results/2018-11-07/{genome}/{sample}_transcripts.txt"
    output:
        "results/2018-11-07/{genome}/{sample}_transcripts.bed"
    conda:
        "../../envs/homer.yaml"
    shell:
        "pos2bed.pl {input} | sort -k1,1 -k2,2n - > {output}"
