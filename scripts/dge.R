## FIXME BiocManager::install("edgeR")

library(readr)
library(edgeR)
library(limma)
library(Homo.sapiens)

## Read in Genes
raw <- read.delim(snakemake@input[[1]],stringsAsFactors=FALSE,check.names=FALSE,)

## Take column 7 and beyond for counts
## First Column is the RefSeqID
y <- DGEList(counts=raw[,7:ncol(raw)], genes=raw[,1])
samplenames <- substring(colnames(y), 12, nchar(colnames(y)))
## FIXME Do dynamically
group <- as.factor(c("Early", "Early", "Early", "Middle", "Middle", "Late", "Late", "Late"))
y$samples$group <- group

cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

## Figure 1
print(snakemake@config[["cell"]])
png(snakemake@output[[1]])
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()

## Calc TMM
y <- calcNormFactors(y, method = "TMM")
y$samples$norm.factors

x2 <- y
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

## Figure 2
png(snakemake@output[[2]])
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data", ylab="Log-cpm")
dev.off()

## Fig 3
## TODO split into groups based on time (early, middle,late)
numColmnSplit <- ncol(raw) / 3
print(numColmnSplit)
## So take total number / 3, and then color each column accordingly
png(snakemake@output[[3]])
lcpm <- cpm(y, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, col=col.group)
title(main="A. Sample groups")
dev.off()

## DGE
## library(splines)
## X <- ns(y$Time, df=5)
## Group <- factor(y$Group)
## design <- model.matrix(~Group*X)
## fit <- lmFit(y, design)
## fit <- eBayes(fit)
## topTable(fit, coef=8:12)

## write_tsv(fit, snakemake@output[[1]])
