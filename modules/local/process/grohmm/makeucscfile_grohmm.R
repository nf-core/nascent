#!usr/bin/env Rscript

#
# REQUIREMENTS
# Packages below need to be available to load when running R.
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID - this will be used for pooling biological replicates.

library(groHMM)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


option_list <- list(
    make_option(c("-o", "--samplesheet"        ), type="character", default='./'    , metavar="path"   , help="Samplesheet with replicates labelled."                                                                      ),
    make_option(c("-i", "--bed_files"    ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bed files."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$bed_files)){
    print_help(opt_parser)
    stop("Please provide bed files and a sample sheet.", call.=FALSE)
}
# Read in sample sheet.
grosamplesheet <- read.delim(file=opt$samplesheet, sep = ",", header=TRUE)
grosamplesheet <- as.data.frame(grosamplesheet)
# Load bed files as well.
# Begin by reading files from sample sheet to id replicates.
S1L0002 <- # add in replicates for merging
# Need to read in strandedness and change to +/-
grosamplesheet$grohmmstrand[grosamplesheet$strandedness == "forward"] <- "+"
grosamplesheet$grohmmstrand[grosamplesheet$strandedness == "reverse"] <- "-"
# Append column with the reverse value
grosamplesheet$grohmmreverse[grosamplesheet$strandedness == "reverse"] <- TRUE
grosamplesheet$grohmmreverse[grosamplesheet$strandedness == "forward"] <- FALSE



if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

# Begin use of groHMM -> NEED TO MAKE REPRODUCIBLE
A <- as(readGAlignments(system.file("extdata", "A.bed", package = "groHMM")), "GRanges")
B <- as(readGAlignments(system.file("extdata", "B.bed", package = "groHMM")), "GRanges")


UCSCfileA <- writeWiggle(reads=opt$outprefix, file= paste(opt$outprefix), "ucscwiggle.wig", fileType="wig", strand="grosamplesheet[1, "strandedness"]"
+ reverse=grosamplesheet[1, grohmmreverse])
\
> # Normalized wiggle files -> NEED TO MAKE REPRODUCIBLE
> expCounts <- mean(c(NROW(A), NROW(B)))

UCSCfileBnormalized <-  writeWiggle(reads=opt$outprefix, file= paste(opt$outprefix), "ucscwigglenormalized.wig", fileType="wig", fileType="wig", strand="+",
+ normCounts=expCounts/NROW(opt$outprefix), reverse=grosamplesheet[1, grohmmreverse])

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
