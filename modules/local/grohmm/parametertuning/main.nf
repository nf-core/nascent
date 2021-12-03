// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '1.24'

process GROHMM_PARAMETERTUNING{
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda    (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-optparse=1.6.6 conda-forge::r-rmariadb=1.1.0 bioconda::bioconductor-edger=3.32.0 bioconda::bioconductor-grohmm=1.24.0 bioconda::bioconductor-org.hs.eg.db=3.12.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2"  : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/mulled-v2-598ad5976af555e8404016e81a380911377ebb95:79c659213212873c7f35071e27b6732b6de63a77-0"
    } else {
    container "quay.io/biocontainers/mulled-v2-598ad5976af555e8404016e81a380911377ebb95:79c659213212873c7f35071e27b6732b6de63a77-0"
    }

    input:
    tuple val(meta), path(bam)
    path(tune)

    output:
    path "*.tuning.csv"    , optional:true    , emit: tuning
    path "*.RData"                    , optional:true    , emit: rdata
    path "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    parameter_tuning.R \\
        --bam_file ${bam} \\
        --tuning_file ${tune} \\
        --outprefix ${prefix} \\
        --genome $params.genome \\
        --outdir ./ \\
        --cores $task.cpus \\
        $options.args

    Rscript -e "library(groHMM); write(x=as.character(packageVersion('groHMM')), file='${software}.version.txt')"
    """
}
