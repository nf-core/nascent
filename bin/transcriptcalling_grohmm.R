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

# Load alignment files
alignments <- c()
for (bam in args$bam_files) {
    alignments <- append(
        alignments,
        as(readGAlignments(bam), "GRanges")
    )
    alignments <- keepStandardChromosomes(alignments, pruning.mode = "coarse")
}

# Call annotations > DEFAULT VALUES ASSIGNED
if (is.null(args$tuning_file)) {
    hmm_result <- detectTranscripts(
        alignments,
        LtProbB = args$ltprobb,
        UTS = args$uts,
        threshold = 1
    ) # Uses either inputted or default values
} else {
    tune <- read.csv(args$tuning_file)
    # Minimum error
    uts <- tune[which.min(tune$errorRate), "UTS"]
    lt_probb <- tune[which.min(tune$errorRate), "LtProbB"]
    hmm_result <- detectTranscripts(
        alignments,
        LtProbB = lt_probb,
        UTS = uts,
        threshold = 1
    )
}

tx_hmm <- hmm_result$transcripts
write.table(
    tx_hmm,
    file = paste(args$outprefix,
        ".transcripts.txt",
        sep = ""
    )
)

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

# Evaluate HMM Annotations
e <- evaluateHMMInAnnotations(tx_hmm, kg_consensus)
# Save as txt file
capture.output(e$eval, file = paste0(args$outprefix, ".eval.txt"))

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
    reads = alignments
)
b_plus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "+")
b_minus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "-")
tx_broken <- c(b_plus, b_minus)
tx_final <- combineTranscripts(tx_broken, kg_consensus)
td_final <- getTxDensity(tx_final, con_expressed, mc.cores = args$cores)
export(
    tx_final,
    con = paste(args$outprefix, ".final.transcripts.bed", sep = "")
)
capture.output(td_final, file = paste0(args$outprefix, ".tdFinal.txt"))
# Output plot
jpeg(file = paste0(args$outprefix, ".tdplot_mqc.jpg"))
# 2. Create the plot
td_final <- getTxDensity(tx_final, con_expressed, mc.cores = args$cores)

# 3. Close the file
dev.off()


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
