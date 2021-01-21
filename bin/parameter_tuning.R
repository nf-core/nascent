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
library(RMariaDB)
library(AnnotationDbi)

option_list <- list(
    make_option(c("-i", "--bam_files"     ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bam files."),
 #   make_option(c("-i", "--ref_transcript"), type="character", default=NULL    , metavar="path"   , help= "Reference transcript annotations."),
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
readsfile <- as(readGAlignments(file = opt$count_file, "GRanges"))

# Call annotations > DEFAULT VALUES ASSIGNED
hmmResult <- detectTranscripts(readsfile, LtProbB=-200, UTS=5, threshold=1)
txHMM <- hmmResult$transcripts
# Input transcript annotations > CURRENTLY JUST USES R LIBRARY > can be changed to generate from UCSC
kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
kgChr7 <- transcripts(kgdb, filter=list(tx_chrom = "chr7"),
columns=c("gene_id", "tx_id", "tx_name"))
seqlevels(kgChr7) <- seqlevelsInUse(kgChr7)
# Collapse annotations in preparation for overlap
kgConsensus <- makeConsensusAnnotations(kgChr7, keytype="gene_id", mc.cores=getOption("mc.cores"))
map <-select(org.Hs.eg.db, keys=unlist(mcols(kgConsensus)$gene_id), columns=c("SYMBOL"), keytype=c("ENTREZID"))
mcols(kgConsensus)$symbol <- map$SYMBOL
mcols(kgConsensus)$type <- "gene"

# TUNING
tune <- data.frame(LtProbB=c(rep(-100,3), rep(-200,3), rep(-300,3)),
UTS=rep(c(5,10,15), 3))
Fp <- windowAnalysis(readsfile, strand="+", windowSize=50)
Fm <- windowAnalysis(readsfile, strand="-", windowSize=50)
evals <- mclapply(seq_len(9), function(x) {
hmm <- detectTranscripts(Fp=Fp, Fm=Fm, LtProbB=tune$LtProbB[x],
UTS=tune$UTS[x])
e <- evaluateHMMInAnnotations(hmm$transcripts, kgConsensus)
e$eval
}, mc.cores=getOption("mc.cores"), mc.silent=TRUE)
tune <- cbind(tune, do.call(rbind, evals))
write.table(tune, file = paste(opt$outprefix,"tuning.parameters.txt", sep=""))


# CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")
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
