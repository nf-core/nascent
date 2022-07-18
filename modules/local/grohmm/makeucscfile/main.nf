process GROHMM_MAKEUCSCFILE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda    (params.enable_conda ? "conda-forge::r-base=4.1.1 conda-forge::r-optparse=1.7.1 conda-forge::r-argparse=2.1.3 bioconda::bioconductor-genomicfeatures=1.46.1 bioconda::bioconductor-grohmm=1.28.0"  : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' :
        'quay.io/biocontainers/mulled-v2-e9a6cb7894dd2753aff7d9446ea95c962cce8c46:0a46dae3241b1c4f02e46468f5d54eadcf64beca-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    path "*.collapsed.wig", emit: wig
    path "*.plus.wig"     , emit: pluswig
    path "*.minus.wig"    , emit: minuswig
    path  "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeucscfile.R \\
        --bam_file ${bam} \\
        --outprefix ${prefix} \\
        --outdir ./ \\
        --cores $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-grohmm: \$(Rscript -e "library(groHMM); cat(as.character(packageVersion('groHMM')))")
    END_VERSIONS
    """
}
