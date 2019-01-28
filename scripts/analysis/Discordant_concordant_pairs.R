
setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)

# Discordant/concordant enhancer list

# Load RPKM matrix
load("Data/eRNA_expression_GM_unfiltered.RData")
gene.rpkm <- read.table("Data/Expression_matrix.symbol.csv")

rpkm <- rbind(eRNA.rpkm, gene.rpkm)

#################################################################################################################

# Heatmaps
plot_heatmap <- function(l_dc, l_cc){
  gene_list <- c(l_dc, l_cc)
  heat.data <- rpkm[gene_list,]
  heat.data <- (heat.data - rowMeans(heat.data)) / apply(heat.data, 1, sd)
  heat.data <- as.matrix(heat.data)
  for(i in c(1:length(gene_list))){
    rownames(heat.data)[i] <- substr(gene_list[i], 14, 1000L)
  }
  for(i in c(1:ncol(heat.data))){
    colnames(heat.data)[i] <- substr(colnames(heat.data)[i], 3, 1000L)
  }
  my_palette <- colorRampPalette(c("blue", "white", "yellow"), bias=1)(n = 100)
  breaks <- c(seq(-3, 3, length.out=101))
  heatmap.2(heat.data, col=my_palette,
            Colv=F, Rowv=F, dendrogram="none", labRow='',
            RowSideColors = c(rep('darkgrey', length(l_dc)), rep('lightgrey', length(l_cc))),
            keysize=1.2, density.info="density", key.ylab="", key.xlab="", key.title="",
            trace="none")
}

dc <- read.table("~/Documents/eRNA/link/IFN project/discordant_pairs_gene-control-enhancer.pairs.txt", header = F, colClasses = 'character')
cc <- read.table("~/Documents/eRNA/link/IFN project/concordant_pairs_gene-control-enhancer.pairs.txt", header = F, colClasses = 'character')
pdf("Figure/Discordant_concordant_pairs_gene.pdf")
plot_heatmap(dc[,1], cc[,1])
dev.off()
pdf("Figure/Discordant_concordant_pairs_eRNA.pdf")
plot_heatmap(dc[,2], cc[,2])
dev.off()
