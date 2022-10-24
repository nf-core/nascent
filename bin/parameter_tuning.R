#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(groHMM))

parser <- ArgumentParser(description = "Run groHMM on some bam files")

parser$add_argument(
    "-i",
    "--bam_files",
    type = "character",
    nargs = "+",
    metavar = "path",
    help = "GRO SEQ data in bam files.",
    required = TRUE
)
parser$add_argument(
    "-t",
    "--tuning_file",
    type = "character",
    default = NULL,
    metavar = "path",
    help = "File with tuning parameters and error rates."
)
parser$add_argument(
    "-o",
    "--outdir",
    type = "character",
    default = "./",
    metavar = "path",
    help = "Output directory."
)
parser$add_argument(
    "-l",
    "--ltprobb",
    type = "integer",
    default = -200,
    metavar = "integer",
    help = cat(
        "Log-transformed transition probability of switching from transcribed
        state to non-transcribed state"
    )
)
parser$add_argument(
    "-u",
    "--uts",
    type = "integer",
    default = 5,
    metavar = "integer",
    help = cat(
        "Variance of the emission probability for reads in the
        non-transcribed state, respectively."
    )
)
parser$add_argument(
    "-p",
    "--outprefix",
    type = "character",
    default = "grohmm",
    metavar = "string",
    help = "Output prefix."
)
parser$add_argument(
    "-g",
    "--gtf",
    type = "character",
    default = NULL,
    metavar = "string",
    help = "GTF File to create TxDb",
    required = TRUE
)
parser$add_argument(
    "-c",
    "--cores",
    type = "integer",
    default = 1,
    metavar = "integer",
    help = "Number of cores."
)

args <- parser$parse_args()

setwd(args$outdir)

if (is.null(args$bam_files)) {
    print_help(args)
    stop("Please provide a bam file", call. = FALSE)
}
if (is.null(args$tuning_file)) {
    print_help(args)
    stop("Please provide a tuning file", call. = FALSE)
}


# Read in bam file.
if (file.exists(args$outdir) == FALSE) {
    dir.create(args$outdir, recursive = TRUE)
}
setwd(args$outdir)

# CHANGE BASED ON PAIRED OR SINGLE END
alignments <- c()
for (bam in args$bam_files) {
    alignments <- append(
        alignments,
        as(readGAlignments(bam), "GRanges")
    )
    alignments <- keepStandardChromosomes(alignments, pruning.mode = "coarse")
}

print("Input transcript annotations")
kg_db <- makeTxDbFromGFF(args$gtf)
kg_tx <- transcripts(kg_db, columns = c("gene_id", "tx_id", "tx_name"))
print("Collapse annotations in preparation for overlap")
kg_consensus <- makeConsensusAnnotations(
    kg_tx,
    keytype = "gene_id",
    mc.cores = args$cores
)
print("Finished consensus annotations")

# TUNING
tune <- read.csv(args$tuning_file)

evals <- mclapply(seq_len(nrow(tune)), function(x) {
    hmm <- detectTranscripts(
        reads = alignments,
        LtProbB = tune$LtProbB[x], UTS = tune$UTS[x]
    )
    e <- evaluateHMMInAnnotations(hmm$transcripts, kg_consensus)
    e$eval
}, mc.cores = args$cores, mc.silent = TRUE)

tune <- cbind(tune, do.call(rbind, evals))
write.csv(tune, file = paste0(args$outprefix, ".tuning.csv"))


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
