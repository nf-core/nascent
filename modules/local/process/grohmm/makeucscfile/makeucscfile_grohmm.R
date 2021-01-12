#!usr/bin/env Rscript
# TODO
# Allow for user directed normalization and allow for biological replicate cases to be considered.
# REQUIREMENTS
# Packages below need to be available to load when running R.


library(groHMM)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(optparse)
library(GenomicAlignments)

option_list <- list(
    make_option(c("-i", "--bam_files"    ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bam files."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='grohmm', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)
print(opt$bam_files)
if (is.null(opt$bam_files)){
    print_help(opt_parser)
    stop("Please provide bam files", call.=FALSE)

}

# Read in bam file.

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

# Begin use of groHMM -> CURRENTLY ONLY TAKES ONE FILE
readsfile <- as(GenomicAlignments::readGAlignments(file = opt$bam_files, use.names = TRUE))
# Generate wig files
writeWiggle(reads= readsfile, file = paste(opt$outprefix, ".fwd.wig"), strand = "+", reverse = FALSE)
writeWiggle(reads= readsfile, file = paste(opt$outprefix, ".fwd.wig"), strand = "-", reverse = TRUE)
# Generate normalized wig files with the number of reads normalizing -> ALLOW FOR USER INPUT
writeWiggle(reads= readsfile, file = paste(opt$outprefix,".fwd.normalized.wig"), strand = "+", reverse = FALSE, norm = NROWS(readsfile))
writeWiggle(reads= readsfile, file = paste(opt$outprefix,".rev.normalized.wig"), strand = "-", reverse = TRUE,  norm = NROWS(readsfile))

## CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")

## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}

################################################################################################
