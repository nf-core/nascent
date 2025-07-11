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
  "--gxf",
  type = "character",
  default = NULL,
  metavar = "string",
  help = "GFF/GTF File to create TxDb",
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
parser$add_argument(
  "-m",
  "--memory",
  type = "integer",
  default = 56000,
  metavar = "integer",
  help = "Amount of memory in MB"
)

args <- parser$parse_args()

options(mc.cores = getCores(args$cores))
## Windows Specific doesn't actually do anything
memory.limit(size = args$memory)
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
  # Use user supplied values or defaults
  hmm_result <- detectTranscripts(
    alignments,
    LtProbB = args$ltprobb,
    UTS = args$uts,
    threshold = 1
  )
} else {
  # Use tuning file and take the minimum values
  tune <- read.csv(args$tuning_file)
  # Minimum error
  uts <- tune[which.min(tune$errorRate), "UTS"]
  lt_probb <- tune[which.min(tune$errorRate), "LtProbB"]
  # Print the minimums for debugging
  cat("minimum uts:", uts)
  cat("minimum lt_probb:", lt_probb)
  hmm_result <- detectTranscripts(
    alignments,
    LtProbB = lt_probb,
    UTS = uts,
    threshold = 1
  )
}

# NOTE TUNING IN A DIFFERENT SCRIPT

tx_hmm <- hmm_result$transcripts
write.table(
  tx_hmm,
  file = paste(args$outprefix,
    ".transcripts.txt",
    sep = ""
  )
)

print("Input transcript annotations")
kg_db <- makeTxDbFromGFF(args$gxf)
kg_tx <- transcripts(kg_db, columns = c("gene_id", "tx_id", "tx_name"))
# TODO I wonder if I could speed things up by filtering
# by chromosome at the Nextflow level...
# https://github.com/google/deepvariant/issues/744
#                         filter=list(tx_chrom="chr7"))
# exclude any transcripts that are located on chromosomes labeled with "random".
kg_tx <- kg_tx[grep("random", as.character(seqnames(kg_tx)), invert = TRUE), ]
print("Printing kg_tx.......")
print(kg_tx)
print("Collapse annotations in preparation for overlap")
kg_consensus <- makeConsensusAnnotations(
  kg_tx,
  mc.cores = args$cores
)
print("Finished consensus annotations")

# Evaluate HMM Annotations
e <- evaluateHMMInAnnotations(tx_hmm, kg_consensus)
# Save as txt file
capture.output(e$eval, file = paste0(args$outprefix, ".eval.txt"))


print("repairing with annotations")
get_expressed_annotations <- function(features, reads) {
  f_limit <- limitToXkb(features)
  count <- countOverlaps(f_limit, reads)
  features <- features[count != 0, ]
  return(features[
    (quantile(width(features), .05) < width(features)) &
      (width(features) < quantile(width(features), .95)),
  ])
}
con_expressed <- get_expressed_annotations(
  features = kg_consensus,
  reads = alignments
)
b_plus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "+")
b_minus <- breakTranscriptsOnGenes(tx_hmm, kg_consensus, strand = "-")
tx_broken <- c(b_plus, b_minus)
# Assign unique IDs if they're missing
if (
  is.null(mcols(tx_broken)$transcript_id) ||
    any(is.na(mcols(tx_broken)$transcript_id))
) {
  mcols(tx_broken)$transcript_id <- paste0("TX", seq_along(tx_broken))
}

# Filter out any transcripts with NA values in start or end positions
tx_broken_filtered <-
  tx_broken[!is.na(start(tx_broken)) & !is.na(end(tx_broken))]

# Ensure that kg_consensus also doesn't contain NA values
kg_consensus_filtered <-
  kg_consensus[!is.na(start(kg_consensus)) & !is.na(end(kg_consensus))]

# Now call combineTranscripts with the filtered data
tx_final <- combineTranscripts(tx_broken_filtered, kg_consensus_filtered)
export(
  tx_final,
  con = paste(args$outprefix, ".final.transcripts.bed", sep = "")
)
# 1. Output plot
png(
  file = paste0(args$outprefix, ".tdplot_mqc.png"),
  width = 800,
  height = 600, res = 300
)

# 2. Create the plot and capture data
td_final <- getTxDensity(tx_final, con_expressed, mc.cores = args$cores)
# 3. Close the file
dev.off()
capture.output(td_final, file = paste0(args$outprefix, ".tdFinal.txt"))

# Write the data used in the plot to a CSV file
data_to_write <- data.frame(x = td_final$x, profile = td_final$profile)
write.csv(data_to_write,
  file = paste0(args$outprefix, ".tdFinal_mqc.csv"),
  row.names = FALSE
)

########################
## CITE PACKAGES USED ##
########################
citation("groHMM")
citation("GenomicFeatures")
citation("GenomicAlignments")
citation("AnnotationDbi")

####################
## R SESSION INFO ##
####################
r_log_file <- "R_sessionInfo.log"
if (file.exists(r_log_file) == FALSE) {
  sink(r_log_file)
  a <- sessionInfo()
  print(a)
  sink()
}
