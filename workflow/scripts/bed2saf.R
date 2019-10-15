saf <- read.table(snakemake@input[[1]])
colnames(saf) <- c("Chr", "Start", "End", "GeneID", "frame", "Strand")
saf <- saf[, c("GeneID", "Chr", "Start", "End", "Strand")]
write.table(saf, snakemake@output[[1]], sep = "\t", col.names = TRUE, row.names = FALSE)
