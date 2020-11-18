#!usr/bin/env Rscript
# TODO
# Allow for user to input their own annotation dataset
# REQUIREMENTS
# Packages below need to be available to load when running R.


library(groHMM)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tigerstats)
library(RMariaDB)
library(AnnotationDbi)

option_list <- list(
    make_option(c("-i", "--bam_files"     ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bam files."),
    make_option(c("-i", "--ref_transcript"), type="character", default=NULL    , metavar="path"   , help= "Reference transcript annotations."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='grohmm', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

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
readsfile <- as(readGAlignments(file = opt$count_file, header = TRUE), "GRanges")

# Call annotations > DEFAULT VALUES ASSIGNED
hmmResult <- detectTranscripts(readsfile, LtProbB=-200, UTS=5, threshold=1)
transcriptcall <- hmmResult$transcripts
write.table(transcriptcall, file = paste(opt$outprefix,".transcripts.txt", sep=""))

# Input transcript annotations > CURRENTLY JUST USES R LIBRARY > can be changed to generate from USCS
refseq <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")
rgseq <- transcripts(rgdb, columns=c("gene_id", "tx_id", "tx_name"))
seqlevels(rgseq) <- seqlevelsInUse(rgseq)

rgConsensus <- makeConsensusAnnotations(rgseq, keytype="gene_id", mc.cores=opt$cores)
# ensure select isn't masked by another package
map <- AnnotationDbi::select(org.Hs.eg.db, keys=unlist(mcols(kgConsensus)$gene_id), columns=c("SYMBOL"), keytype=c("ENTREZID"))
mcols(rgConsensus)$symbol <- map$SYMBOL
mcols(rgConsensus)$type <- "gene"

# Evaluate HMM Annotations
e <- evaluateHMMInAnnotations(transcriptcall, rgConsensus)
# Save as txt file
capture.output(e$eval, file = paste(opt$outprefix, ".eval.txt", header = TRUE))

## CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")
citation("tigerstats")
citation("RMariaDB")
citation("AnnotationDbi")

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
