#!/usr/bin/env Rscript
# TODO
# Allow for user to input their own annotation dataset
# REQUIREMENTS
# Packages below need to be available to load when running R.


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(groHMM))

option_list <- list(
    make_option(c("-i", "--bam_file"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "Time course of GRO SEQ data in bam files."
    ),
    make_option(c("-t", "--tuning_file"),
        type = "character",
        default = NULL,
        metavar = "path",
        help = "File with all possible LTS/UTS combinations."
    ),
    make_option(c("-o", "--outdir"),
        type = "character",
        default = "./",
        metavar = "path",
        help = "Output directory."
    ),
    make_option(c("-g", "--genome"),
        type = "character",
        default = "hg19",
        metavar = "string",
        help = "Reference UCSC genome"
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
if (is.null(opt$tuning_file)) {
    print_help(opt_parser)
    stop("Please provide a tuning file", call. = FALSE)
}


# Read in bam file.
if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive = TRUE)
}
setwd(opt$outdir)

# CHANGE BASED ON PAIRED OR SINGLE END
galigned <- readGAlignments(BamFile(opt$bam_file, asMates = TRUE))
reads_file <- GRanges(galigned)

# Call annotations > DEFAULT VALUES ASSIGNED
hmm_result <- detectTranscripts(
    reads_file,
    LtProbB = -200,
    UTS = 5,
    threshold = 1
)
tx_hmm <- hmm_result$transcripts
write.table(tx_hmm, file = paste(opt$outprefix, ".transcripts.txt", sep = ""))

kg_db <- makeTxDbFromUCSC(genome = opt$genome, tablename = "refGene")
saveDb(kg_db, file = "RefGene.sqlite")
kg_db <- loadDb("RefGene.sqlite")
kg_tx <- transcripts(kg_db, columns = c("gene_id", "tx_id", "tx_name"))
# Collapse annotations in preparation for overlap
kg_consensus <- makeConsensusAnnotations(
    kg_tx,
    keytype = "gene_id",
    mc.cores = opt$cores
)

map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unlist(mcols(kg_consensus)$gene_id),
    columns = c("SYMBOL"),
    keytype = c("ENTREZID")
)
mcols(kg_consensus)$symbol <- map$SYMBOL
mcols(kg_consensus)$type <- "gene"

# TUNING
tune <- read.csv(opt$tuning_file)

# Default values
## tune <- data.frame(
##   LtProbB = c(
##     rep(-100, 9),
##     rep(-150, 9),
##     rep(-200, 9),
##     rep(-250, 9),
##     rep(-300, 9),
##     rep(-350, 9),
##     rep(-400, 9)
##   ),
##   UTS = rep(
##     c(5, 10, 15, 20, 25, 30, 35, 40, 45),
##     7
##   )
## )

evals <- mclapply(seq_len(nrow(tune)), function(x) {
    hmm <- detectTranscripts(
        reads = reads_file,
        LtProbB = tune$LtProbB[x], UTS = tune$UTS[x]
    )
    e <- evaluateHMMInAnnotations(hmm$transcripts, kg_consensus)
    e$eval
}, mc.cores = opt$cores, mc.silent = TRUE)

tune <- cbind(tune, do.call(rbind, evals))
write.csv(tune, file = paste0(opt$outprefix, ".tuning.csv"))


# CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("GenomicAlignments")
citation("AnnotationDbi")

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
