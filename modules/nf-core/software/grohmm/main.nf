// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '1.24'

process GROHMM_MAKEUCSCFILE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-optparse=1.6.6 bioconda::bioconductor-genomic-features=1.32.3 bioconda::bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2 bioconda::bioconductor-edgeR=3.32.0 bioconda::bioconductor-org.hs.eg.db=3.6.0 bioconda::bioconductor-grohmm=1.24.0" : null)
    // TODO build container using mulled toolkit
    container "quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"


    input:
    tuple val(meta), path(bam)

    output:
    path "*fwd.wig"                  , optional:true    , emit: fwdwig
    path "*rev.wig"                  , optional:true    , emit: revwig
    path "*fwd.normalized.wig"       , optional:true    , emit: fwdwig_normalized
    path "*rev.normalized.wig"       , optional:true    , emit: revwig_normalized
    path "*.RData"                   , optional:true    , emit: rdata
    path "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    makeucscfile_grohmm.R \\

        --bam_files $bam \\
        --outdir ./ \\
        --cores $task.cpus \\
        $options.args

    Rscript -e "library(groHMM); write(x=as.character(packageVersion('groHMM')), file='${software}.version.txt')"

    """
}
