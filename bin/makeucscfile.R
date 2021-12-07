#!/usr/bin/env Rscript
# TODO
# Allow for biological replicate cases to be considered.
# REQUIREMENTS
# Packages below need to be available to load when running R.


suppressPackageStartupMessages(library(groHMM))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(RMariaDB))

option_list <- list(
    make_option(c("-i", "--bam_file"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "Time course of GRO SEQ data in bam files."
    ),
    make_option(c("-o", "--outdir"),
        type = "character",
        default = "./",
        metavar = "path",
        help = "Output directory."
    ),
    make_option(c("-n", "--norm"),
        type = "integer",
        default = 1,
        metavar = "integer",
        help = "A normalization factor correcting for library size or other effects. For example, total mappible read counts might be a reasonable value. Default: 1 (i.e. no normalization)."
    ),
    make_option(c("-p", "--outprefix"),
        type = "character",
        default = "grohmm",
        metavar = "string",
        help = "Output prefix."
    ),
    make_option(c("-c", "--cores"),
        type = "integer",
        default = 1,
        metavar = "integer",
        help = "Number of cores."
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt$bam_file)

if (is.null(opt$bam_file)) {
    print_help(opt_parser)
    stop("Please provide a bam file", call. = FALSE)
}

# Read in bam file.


if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive = TRUE)
}
setwd(opt$outdir)

# Begin use of groHMM
# TODO CURRENTLY ONLY TAKES ONE FILE
# TODO CHANGE BASED ON PAIRED OR SINGLE END
galigned <- readGAlignments(BamFile(opt$bam_file, asMates = TRUE))
readsfile <- GRanges(galigned)

# Generate normalized wig files with the number of reads normalizing
# TODO ALLOW FOR USER INPUT, default has no normalization
writeWiggle(
    reads = readsfile,
    file = paste0(
        opt$outprefix,
        ".plus.wig"
    ),
    strand = "+",
    norm = opt$norm
)
writeWiggle(
    reads = readsfile,
    file = paste0(opt$outprefix, ".minus.wig"),
    strand = "-",
    norm = opt$norm
)
writeWiggle(
    reads = readsfile,
    file = paste0(opt$outprefix, ".collapsed.wig"),
    strand = "*",
    norm = opt$norm
)

## CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")

## R SESSION INFO                             ##
################################################
################################################

r_log_file <- "R_sessionInfo.log"
if (file.exists(r_log_file) == FALSE) {
    sink(r_log_file)
    a <- sessionInfo()
    print(a)
    sink()
}

################################################################################
