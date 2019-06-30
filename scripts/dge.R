## FIXME BiocManager::install("edgeR")

library(readr)
library(edgeR)
library(limma)
library(org.Hs.eg.db)

## Read in Genes
raw <- read.delim(snakemake@input[[1]],stringsAsFactors=FALSE,check.names=FALSE,)

## Take column 7 and beyond for counts
## First Column is the RefSeqID
y <- DGEList(counts=raw[,7:ncol(raw)], genes=raw[,1])
head(y$genes)

## Keep only transcripts with IDs in NCBI
## idfound <- y$genes %in% mappedRkeys(org.Hs.egREFSEQ)
## y <- y[idfound,]
## head(y$genes)

print(0)
egREFSEQ <- toTable(org.Hs.egREFSEQ)
print(2)
m <- match(y$genes, egREFSEQ$accession)
print(3)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
print(4)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

## Remove low counts
print(5)
keep <- filterByExpr(y, design)
print(6)
y <- y[keep,,keep.lib.sizes=FALSE]
print(7)
dge <- calcNormFactors(dge)

## limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

write_tsv(fit, snakemake@output[[1]])
