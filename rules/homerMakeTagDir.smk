rule All_makeTagDirectory:
    input:
        expand("results/2018-10-04/{unit}.bam",unit=SAMPLES)
    output:
        "results/2018-11-07/All_together/"
    conda:
        "../envs/homer.yaml"
    shell:
        "makeTagDirectory {output} -genome hg19 -checkGC {input}"
