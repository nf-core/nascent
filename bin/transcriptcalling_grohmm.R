#!/usr/bin/env Rscript
# TODO
# Allow for user to input their own annotation dataset
# REQUIREMENTS
# Packages below need to be available to load when running R.


library(groHMM)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(edgeR)
library(optparse)
library(GenomicAlignments)
library(RMariaDB)

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
        help = "File with tuning parameters and error rates."
    ),
    make_option(c("-o", "--outdir"),
        type = "character",
        default = "./",
        metavar = "path",
        help = "Output directory."
    ),
    make_option(c("-l", "--ltprobb"),
        type = "integer",
        default = -200,
        metavar = "integer",
        help = "Log-transformed transition probability of switching from transcribed state to non-transcribed state"
    ),
    make_option(c("-u", "--uts"),
        type = "integer",
        default = 5,
        metavar = "integer",
        help = "Variance of the emission probability for reads in the non-transcribed state, respectively."
    ),
    make_option(c("-p", "--outprefix"),
        type = "character",
        default = "grohmm",
        metavar = "string",
        help = "Output prefix."
    ),
    make_option(c("-g", "--gtf"),
        type = "character",
        default = NULL,
        metavar = "string",
        help = "GTF File to create TxDb"
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
## readsfile <- as(
##     GenomicAlignments::readGAlignments(
##         file = opt$bam_file,
##         use.names = TRUE
##     )
## )
# TODO CHANGE BASED ON PAIRED OR SINGLE END
galigned <- readGAlignments(BamFile(opt$bam_file, asMates = TRUE))
reads_file <- GRanges(galigned)

# Call annotations > DEFAULT VALUES ASSIGNED
if (is.null(opt$tuning_file)) {
    hmm_result <- detectTranscripts(
        reads_file,
        LtProbB = opt$ltprobb,
        UTS = opt$uts,
        threshold = 1
    ) # Uses either inputted or default values
} else {
    tune <- read.csv(opt$tuning_file)
    # Minimum error
    uts <- tune[which.min(tune$total), "UTS"]
    lt_probb <- tune[which.min(tune$total), "LtProbB"]
    hmm_result <- detectTranscripts(
        reads_file,
        LtProbB = lt_probb,
        UTS = uts,
        threshold = 1
    )
}

tx_hmm <- hmm_result$transcripts
write.table(
    tx_hmm,
    file = paste(opt$outprefix,
        ".transcripts.txt",
        sep = ""
    )
)

# Input transcript annotations
kg_db <- makeTxDbFromGFF(opt$gtf)
kg_tx <- transcripts(kg_db, columns = c("gene_id", "tx_id", "tx_name"))
# Collapse annotations in preparation for overlap
kg_consensus <- makeConsensusAnnotations(
    kg_tx,
    keytype = "gene_id",
    mc.cores = opt$cores
)

# Evaluate HMM Annotations
e <- evaluateHMMInAnnotations(tx_hmm, kg_consensus)
# Save as txt file
capture.output(e$eval, file = paste0(opt$outprefix, ".eval.txt"))

# TUNING IN A DIFFERENT SCRIPT

# repairing with annotations
get_expressed_annotations <- function(features, reads) {
    f_limit <- limitToXkb(features)
    count <- countOverlaps(f_limit, reads)
    features <- features[count != 0, ]
    return(features[(quantile(width(features), .05) < width(features)) &
        (width(features) < quantile(width(features), .95)), ])
}
con_expressed <- get_expressed_annotations(
    features = kg_consensus,
    reads = reads_file
)
b_plus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "+")
b_minus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "-")
tx_broken <- c(b_plus, b_minus)
tx_final <- combineTranscripts(tx_broken, kg_consensus)
td_final <- getTxDensity(tx_final, con_expressed, mc.cores = opt$cores)
export(tx_final, con = paste(opt$outprefix, "final.transcripts.bed", sep = ""))
capture.output(td_final, file = paste0(opt$outprefix, ".tdFinal.txt"))
# Output plot
jpeg(file = paste0(opt$outprefix, ".tdplot.jpg"))
# 2. Create the plot
td_final <- getTxDensity(tx_final, con_expressed, mc.cores = opt$cores)

# 3. Close the file
dev.off()


# CITE PACKAGES USED
citation("groHMM")
citation("GenomicFeatures")
citation("org.Hs.eg.db")
citation("edgeR")
citation("TxDb.Hsapiens.UCSC.hg19.knownGene")
# citation("RMariaDB")
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
