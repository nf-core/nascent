#!/usr/bin/env Rscript
# TODO
# Allow for user to input their own annotation dataset
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
    make_option(c("-i", "--bam_file"      ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bam files."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-l", "--ltprobb" ), type="integer", default=-200         , metavar="integer", help="Log-transformed transition probability of switching from transcribed state to non-transcribed state"                                                                  ),
    make_option(c("-u", "--uts"         ), type="integer", default=5         , metavar="integer", help="Variance of the emission probability for reads in the non-transcribed state, respectively."                                                                  ),
    make_option(c("-p", "--outprefix"     ), type="character", default='grohmm', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

print(opt$bam_file)

if (is.null(opt$bam_file)){
    print_help(opt_parser)
    stop("Please provide a bam file", call.=FALSE)

}

# Read in bam file.


if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

# Begin use of groHMM -> CURRENTLY ONLY TAKES ONE FILE
#readsfile <- as(GenomicAlignments::readGAlignments(file = opt$bam_file, use.names = TRUE))
galigned <- readGAlignments(BamFile(opt$bam_file, asMates=TRUE)) # CHANGE BASED ON PAIRED OR SINGLE END
readsfile <- GRanges(galigned)

#  make_option(c("-i", "--ref_transcript"), type="character", default=NULL    , metavar="path"   , help= "Reference transcript annotations."),

# Call annotations > DEFAULT VALUES ASSIGNED
hmmResult <- detectTranscripts(readsfile, LtProbB=opt$ltprobb, UTS=opt$uts, threshold=1)
txHMM <- hmmResult$transcripts
write.table(txHMM, file = paste(opt$outprefix,".transcripts.txt", sep=""))
# TODO make reproducible, ask for sample file
# Input transcript annotations > CURRENTLY JUST USES R LIBRARY > can be changed to generate from UCSC
kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
kgtx <- transcripts(kgdb, columns=c("gene_id", "tx_id", "tx_name"))
# Collapse annotations in preparation for overlap
kgConsensus <- makeConsensusAnnotations(kgtx, keytype="gene_id", mc.cores=opt$cores)
map <-select(org.Hs.eg.db, keys=unlist(mcols(kgConsensus)$gene_id), columns=c("SYMBOL"), keytype=c("ENTREZID"))
mcols(kgConsensus)$symbol <- map$SYMBOL
mcols(kgConsensus)$type <- "gene"

# Evaluate HMM Annotations
e <- evaluateHMMInAnnotations(txHMM, kgConsensus)
# Save as txt file
capture.output(e$eval, file = paste(opt$outprefix, ".eval.txt", header = TRUE))

# TUNING IN A DIFFERENT SCRIPT

#repairing with annotations
getExpressedAnnotations <- function(features, reads) {
fLimit <- limitToXkb(features)
count <- countOverlaps(fLimit, reads)
features <- features[count!=0,]
return(features[(quantile(width(features), .05) < width(features))
& (width(features) < quantile(width(features), .95)),])}
conExpressed <- getExpressedAnnotations(features=kgConsensus,reads=readsfile)
bPlus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="+")
bMinus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="-")
txBroken <- c(bPlus, bMinus)
txFinal <- combineTranscripts(txBroken, kgConsensus)
tdFinal <- getTxDensity(txFinal, conExpressed, mc.cores=opt$cores)
write.table(txFinal, file = paste(opt$outprefix,"final.transcripts.txt", sep=""))
capture.output(tdFinal, file = paste(opt$outprefix, ".tdFinal.txt", header = TRUE))
#Output plot
jpeg(file = paste(opt$outprefix, ".tdplot.jpg", header = TRUE))
# 2. Create the plot
tdFinal <- getTxDensity(txFinal, conExpressed, mc.cores=opt$cores)
# 3. Close the file
dev.off()


# CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")
#citation("RMariaDB")
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
