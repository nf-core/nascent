library(edgeR)
library(org.Hs.eg.db)

fc <- readRDS(file = snakemake@input[[1]])

group <- colnames(fc$counts)
y <- DGEList(counts=fc$counts, group=group, genes=fc$annotation)

Symbol <- mapIds(org.Hs.eg.db, keys=rownames(y), keytype="ENTREZID", column="SYMBOL")
y$genes <- data.frame(Symbol=Symbol)

## Tae's
rawdata_table <- y$counts
genelen <- cbind(y$genes$GeneID, y$genes$Length)
write.csv(rawdata_table, file=snakemake@output[[1]]) #you can use this file in PIVOT
write.csv(genelen, file=snakemake@output[[2]]) #you need this file in PIVOT for normalization

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
## write.csv(y$counts, file=snakemake@output[[3]]) #you can use this file in PIVOT


bcv <- 0.2
et <- exactTest(y, dispersion=bcv^2)
## head(et)
topGenes <- topTags(et)
print(topGenes)

## Common dispersion
## y1 <- y
## y1$samples$group <- 1
## y0 <- estimateDisp(y1[topGenes,], trend="none", tagwise=FALSE)
## y$common.dispersion <- y0$common.dispersion
## fit <- glmFit(y, design)
## lrt <- glmLRT(fit)
