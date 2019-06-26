library(tidyr)
library(edgeR)
library(limma)
library(org.Hs.eg.db)

## Read in Genes
raw <- read_tsv(snakemake@input[[1]])

## Take column 7 and beyond for counts
## First Column is the RefSeqID
y <- DGEList(counts=raw[,7:ncol(raw)], genes=raw[,1])
head(y)

## Keep only transcripts with IDs in NCBI
idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y <- y[idfound,]

egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
egSYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

## Remove low counts
keep <- filterByExpr(y, design)
y <- y[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

## limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

write_tsv(fit, snakemake@output[[1]])
