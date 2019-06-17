library(biomaRt)
library(dplyr)

counts <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)

human_mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("external_gene_name", "percentage_gene_gc_content", "chromosome_name", "start_position", "end_position", "transcript_biotype"),
    filters = "external_gene_name", values = counts$Geneid, mart = human_mart)
genemap <- distinct(genemap, external_gene_name, .keep_all = TRUE )
rownames(genemap) <- genemap$external_gene_name
genemap$external_gene_name <- NULL
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
mylength <- select(counts, Length)
mygc <- select(genemap, percentage_gene_gc_content)
mybiotypes <- select(genemap, transcript_biotype)
mychroms <- select(genemap, chromosome_name, start_position, end_position)
#mychroms <- select(counts, chrom, txstart, txend)
sampleInfo <- data.frame(
  sample = c("GM0h", "GM12h", "GM18h", "GM1h", "GM24h", "GM2h", "GM30min", "GM48h", "GM4h", "GM6h", "GM72h", "GM9h"),
  treatment = c("0h", "12h", "18h", "1h", "24h", "2h", "30min", "48h", "4h", "6h", "72h", "9h"))
library(NOISeq)
mycounts <- counts[, 6:17]
mydata <- readData(data = mycounts, length = mylength, chromosome = mychroms, factors = sampleInfo)
myresults <- noiseq(mydata, k=0.5, norm = "tmm", factor = "treatment", conditions = c("0h", "12h"),  pnr = 0.2, nss = 5, v = 0.02, lc = 0, replicates = "no")
#iterate through all pairs
deg <- degenes(myresults, q = 0.9, M = NULL) #iterate for all results objects
write.table(as.data.frame(deg), file = snakemake@output[[1]], sep = " ", row.names = TRUE, col.names = TRUE)
