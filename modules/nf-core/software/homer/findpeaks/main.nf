// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '4.11'

process HOMER_FINDPEAKS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3 bioconda::samtools=1.11 conda-forge::r-base=4.0.2 bioconda::bioconductor-deseq2=1.30.0 bioconda::bioconductor-edger=3.32.0 conda-forge::perl=5.26.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0"
    }

    input:
    tuple val(meta), path(tagDir)

    output:
    tuple val(meta), path("*peaks.txt"), emit: txt
    path  "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    findPeaks \\
        $tagDir \\
        $options.args \\
        -o ${prefix}.peaks.txt

    echo $VERSION > ${software}.version.txt
    """
}
