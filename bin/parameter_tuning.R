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
memory.limit(size = args$memory)
setwd(args$outdir)

if (is.null(args$bam_files)) {
  print_help(args)
  stop("Please provide a bam file", call. = FALSE)
}

# Load alignment files
# TODO? CHANGE BASED ON PAIRED OR SINGLE END
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
  mc.cores = args$cores
)
print("Finished consensus annotations")

############
## TUNING ##
############
print("Starting tuning run")
Fp <- windowAnalysis(alignments, strand = "+", windowSize = 50)
Fm <- windowAnalysis(alignments, strand = "-", windowSize = 50)
hmm <- detectTranscripts(
  Fp = Fp,
  Fm = Fm,
  reads = alignments,
  LtProbB = args$ltprobb,
  UTS = args$uts
)
print("Evaluating")
e <- evaluateHMMInAnnotations(hmm$transcripts, kg_consensus)
print(e$eval)
print(e)

write.csv(e$eval, file = paste0(args$outprefix, ".tuning.csv"))

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
