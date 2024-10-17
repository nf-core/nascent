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
# Import the GTF file using rtracklayer
gtf <- import(args$gtf)

# Exclude any transcripts located on chromosomes labeled with "random"
gtf <- gtf[!grepl("random", seqnames(gtf)), ]

# Extract transcript-level features
transcripts_gtf <- gtf[gtf$type == "transcript", ]
# Extract exon features
exons_gtf <- gtf[gtf$type == "exon", ]

# Ensure that the 'transcript_id' and 'gene_id' columns are present
if (!all(c("transcript_id", "gene_id") %in% names(mcols(exons_gtf)))) {
  stop("The GTF file lacks 'transcript_id' or 'gene_id' in its attributes.")
}

# Group exons by transcript_id
exons_by_transcript <- split(exons_gtf, exons_gtf$transcript_id)

# Diagnostic prints
print(paste("Number of transcripts:", length(exons_by_transcript)))

# Reduce exons to create transcript ranges
transcripts_ranges <- GenomicRanges::reduce(exons_by_transcript)
transcripts_ranges <- unlist(transcripts_ranges, use.names = TRUE)

# Diagnostic prints after reduction
print(paste("Number of transcripts_ranges after reduction:", length(transcripts_ranges)))

# Create mapping dataframe
mapping_df <- data.frame(
  transcript_id = names(transcripts_ranges),
  gene_id = vapply(exons_by_transcript[names(transcripts_ranges)], function(x) unique(x$gene_id)[1], character(1)),
  stringsAsFactors = FALSE
)

# Check for length mismatch
if (nrow(mapping_df) != length(transcripts_ranges)) {
  stop(paste("Length mismatch between mapping_df and transcripts_ranges:", nrow(mapping_df), length(transcripts_ranges)))
}

# Assign metadata
mcols(transcripts_ranges)$transcript_id <- mapping_df$transcript_id
mcols(transcripts_ranges)$gene_id <- mapping_df$gene_id

# Assign seqnames and strand from the exons
seqnames(transcripts_ranges) <- seqnames(exons_gtf[match(names(transcripts_ranges), exons_gtf$transcript_id)])
strand(transcripts_ranges) <- strand(exons_gtf[match(names(transcripts_ranges), exons_gtf$transcript_id)])

# Ensure that seqlevels are set correctly
seqlevels(transcripts_ranges) <- seqlevels(gtf)

# Remove any transcripts with NA values
transcripts_ranges <- transcripts_ranges[!is.na(start(transcripts_ranges)) & !is.na(end(transcripts_ranges))]

print("Collapse annotations in preparation for overlap")

# Filter out non-standard chromosomes
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
transcripts_ranges_filtered <- transcripts_ranges[seqnames(transcripts_ranges) %in% standard_chromosomes]

# Sort the ranges (important for efficient processing)
transcripts_ranges_filtered <- sort(transcripts_ranges_filtered)

# Process in chunks
chunk_size <- 50000  # Adjust based on your available memory
num_chunks <- ceiling(length(transcripts_ranges_filtered) / chunk_size)

kg_consensus <- GRanges()
for (i in 1:num_chunks) {
  print(paste("Processing chunk", i, "of", num_chunks))
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(transcripts_ranges_filtered))
  chunk <- transcripts_ranges_filtered[start_idx:end_idx]

  kg_consensus_chunk <- makeConsensusAnnotations(
    chunk,
    mc.cores <- min(args$cores, 4) # Adjust based on your HPC setup
  )

  # Preserve gene_id or other relevant metadata
  if ("gene_id" %in% names(mcols(chunk))) {
      mcols(kg_consensus_chunk)$gene_id <- mcols(chunk)$gene_id[1]
  }

  kg_consensus <- c(kg_consensus, kg_consensus_chunk)

  # Force garbage collection to free up memory
  gc()
}

# Final reduction step to ensure non-overlapping annotations
kg_consensus <- reduce(kg_consensus, ignore.strand = FALSE)

# If gene_id was preserved, ensure it's not lost in the reduction
if ("gene_id" %in% names(mcols(kg_consensus))) {
  kg_consensus <- kg_consensus[!is.na(kg_consensus$gene_id)]
}

print("Finished consensus annotations for CHM13")

############
## TUNING ##
############
print("Starting tuning run")
tune <- data.frame(
  LtProbB = args$ltprobb,
  UTS = args$uts
)
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

# Extract evaluation metrics and convert to a data frame
eval_metrics <- as.data.frame(e$eval)

# If eval_metrics is a list of lists, unlist it
if (is.list(eval_metrics[[1]])) {
  eval_metrics <- as.data.frame(t(sapply(e$eval, unlist)))
}

# Combine the tuning parameters with the evaluation metrics
tune <- cbind(tune, eval_metrics)

print(e$eval)
print(e)

# Write the combined data to a CSV file without row names
write.csv(tune, file = paste0(args$outprefix, ".tuning.csv"), row.names = FALSE)
# Write kg_consensus to a bed file for testing
export.bed(kg_consensus, con = paste0(args$outprefix, ".tuning.consensus.bed"))

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
