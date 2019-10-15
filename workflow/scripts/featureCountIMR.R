library(Rsubread)

fc <- featureCounts(snakemake@input[["bams"]], annot.inbuilt="hg19", nthreads=snakemake@threads)

output <- sub(".bam", "", snakemake@input[["bams"]])
outputNames <- sub("results/2018-12-01/IMR/", "", output)

colnames(fc$counts) <- outputNames
head(fc$counts)

saveRDS(fc, file = snakemake@output[[1]])
