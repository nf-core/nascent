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
library(optparse)
library(GenomicAlignments)


option_list <- list(
    make_option(c("-i", "--bam_file"     ), type="character", default=NULL    , metavar="path"   , help="Time course of GRO SEQ data in bam files."),
 #   make_option(c("-i", "--ref_transcript"), type="character", default=NULL    , metavar="path"   , help= "Reference transcript annotations."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='grohmm', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

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
readsfile <- as(readGAlignments(file = opt$count_file, "GRanges"))

galigned <- readGAlignments(BamFile(opt$bam_file, asMates=TRUE)) # CHANGE BASED ON PAIRED OR SINGLE END
readsfile <- GRanges(galigned)

#  make_option(c("-i", "--ref_transcript"), type="character", default=NULL    , metavar="path"   , help= "Reference transcript annotations."),

# Call annotations > DEFAULT VALUES ASSIGNED
hmmResult <- detectTranscripts(readsfile, LtProbB=-200, UTS=5, threshold=1)
txHMM <- hmmResult$transcripts
write.table(txHMM, file = paste(opt$outprefix,".transcripts.txt", sep=""))
# TODO make reproducible, ask for sample file
# Input transcript annotations > CURRENTLY JUST USES R LIBRARY > can be changed to generate from UCSC
kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
kgtx <- transcripts(kgdb, columns=c("gene_id", "tx_id", "tx_name"))
# Collapse annotations in preparation for overlap
kgConsensus <- makeConsensusAnnotations(kgtx, keytype="gene_id", mc.cores=opt$cores)


# TUNING

tune <- data.frame(LtProbB=c(rep(-100,9), rep(-150,9), rep(-200,9), rep(-250,9), rep(-300,9),rep(-350,9),rep(-400,9) ),
           UTS=rep(c(5,10,15,20,25,30,35,40,45), 7))

evals <- mclapply(seq_len(63), function(x) {
        hmm <- detectTranscripts(reads=readsfile, LtProbB=tune$LtProbB[x], UTS=tune$UTS[x])
        e <- evaluateHMMInAnnotations(hmm$transcripts, kgConsensus)
        e$eval
        }, mc.cores=opt$cores,  mc.silent=TRUE)

tune <- cbind(tune, do.call(rbind, evals))
write.table(tune, file = paste(opt$outprefix,".tuning.tsv", sep="\t", row.names=F, col.names=T))


# CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")
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
