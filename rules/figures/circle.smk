rule eRNA_circle:
    input:
        "Input File",
    output:
        "Output File"
    log:
        "logs/figures/eRNA_circle.log"
    threads: 2
    script:
        "../../scripts/figures/circle.py"
