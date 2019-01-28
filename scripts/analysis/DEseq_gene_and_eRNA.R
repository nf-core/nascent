setwd("~/Documents/IFN_enhancer/R/")
library(GenomicRanges)
library(GenomicFeatures)
library(edgeR)
library(DESeq2)
library(cluster)
library(gplots)
source("~/Documents/IFN_enhancer/R/Script/ref_2_symbol_matrix.R")

load("Data/DEseq_input.RData")
#load("eRNA_expression_GM.RData")
#load("gene.rpkm.symbol.RData")

###############################################################
# 1. Perform DEseq
###############################################################
deseq_padj_cutoff <- 0.05
de_times_cutoff <- 1

deg_list <- NULL
deg_ct <- NULL
tp <- unlist(strsplit(colnames(countData), "_"))[c(0:11*4+1)]
groups <- factor(rep(tp, each=2))
colData <- data.frame(groups=groups, row.names = colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ groups)
dds <- dds[ rowSums(rpkmData > 1) > 1, ] # Pre-filtering
dds <- DESeq(dds, betaPrior=FALSE )
norm.counts <- counts(dds, normalized = T)

deg_matrix <- matrix(0, nrow=nrow(countData), ncol=11, dimnames = list(rownames(countData), as.character(tp)[2:12]))
for(i in 2:12){
  res <- results(dds, contrast=c("groups",as.character(tp[i]),"GM0h"))
  lab <- rownames(norm.counts)[which(res$padj < deseq_padj_cutoff)]
  deg_matrix[lab, i-1] <- 1
}

de_all <- rownames(deg_matrix)[rowSums(deg_matrix[,1:9])>de_times_cutoff]
lab <- which(substring(de_all, 1, 12) == "MetaEnhancer")
de_genes <- de_all[-lab]
de_eRNAs <- de_all[lab]

save(norm.counts, countData, deg_matrix, file="Data/deg_matrix_GM.RData")

###############################################################
# 2. Some stats for genes and for eRNAs
###############################################################
load(file="Data/deg_matrix_GM.RData")

deg_matrix <- deg_matrix[rownames(deg_matrix) %in% rownames(norm.counts), ]
deg_matrix_direction <- deg_matrix
ref_col <- c(1:2)
initial <- log2(rowMeans(norm.counts[,ref_col]+1))
for(i in 1:ncol(deg_matrix)){
  lab <- which(deg_matrix[,i]==1)
  a <- log2(rowMeans(norm.counts[lab, i*2 + ref_col])+1)
  b <- initial[lab]
  deg_matrix_direction[lab[a-b<0], i] <- -1
}
count_deg <- function(x){return(c(
  sum(x==1), -sum(x==-1)
))}

save(deg_matrix_direction, file="Data/deg_matrix_direction.RData")

load("Data/deg_matrix_direction.RData")
pdf("Figure/DE_gene_and_eRNA.pdf")
lab <- which(substring(rownames(deg_matrix_direction), 1, 12) == "MetaEnhancer")
# Plot for genes
tp <- apply(deg_matrix_direction[-lab,], 2, count_deg)
barplot(tp[1,], ylim=c(-1200, 1000), col="red", ylab="Number of Differentially Expressed Genes",
        names.arg = substr(colnames(tp), 3, 1000L))
barplot(tp[2,], add=T, col="blue", names.arg = F)

# Plot for enhancers
tp <- apply(deg_matrix_direction[lab,], 2, count_deg)
barplot(tp[1,], ylim=c(-200, 600), col="red", ylab="Number of Differentially Expressed Enhancers",
        names.arg = substr(colnames(tp), 3, 1000L))
barplot(tp[2,], add=T, col="blue", names.arg = F)
dev.off()









###############################################################
# 3. Summarize gene expression profiles near differentially expressed eRNAs
###############################################################

x <- eRNA.rpkm[, 2:10]
lab <- de_eRNAs[de_eRNAs %in% rownames(x)]
x <- x[lab, ]
x <- (x - rowMeans(x)) / apply(x, 1, sd)

# Sort rows by kmeans clustering
zscore <- x
set.seed(0)
centers <- 3
kmeans <- kmeans(zscore, centers, iter.max=10^5)

k_cluster <- kmeans$cluster
cluster_order <- c(2,1,3)
ind <- cluster_order[k_cluster] * 10^4 + rowMeans(zscore)
od <- order(ind, decreasing = T)
zscore <- zscore[od, ]
k_cluster <- k_cluster[od]

# eRNA heatmap
my_palette <- colorRampPalette(c("blue", "black", "yellow"), bias=1)(n = 500)
heatmap.2(zscore[nrow(zscore):1,], col=my_palette, Colv=F, Rowv=F, dendrogram="none",
          labRow=NA, keysize=1.5, density.info="density",
          key.ylab="", key.xlab="", key.title="", trace="none")

par(mfcol=c(2,2), mar=c(4,4,4,5))
for(i in cluster_order){
  boxplot(zscore[k_cluster==i,], ylab="Zscore")
}

###############################################################
# 4. Finding correlated eRNA-gene pairs
###############################################################

tss <- promoters(gene, upstream = 1, downstream = 1)
flanking <- 200000
points_to_consider <- c(1:10)
e_mat <- eRNA.rpkm[, points_to_consider]
g_mat <- gene.rpkm.symbol[, points_to_consider]

find_target <- function(x){
  tp <- eRNA[mcols(eRNA)$id == x]
  start(tp) <- start(tp) - flanking
  end(tp) <- end(tp) + flanking
  hits <- findOverlaps(tp, tss, ignore.strand=T)
  lab <- as.character(unique(mcols(tss[subjectHits(hits)])$id))
  lab <- ref_2_symbol(lab)
  lab <- lab[lab %in% de_genes]
  return(lab)
}

find_region <- function(x, y){
  a <- x
  b <- y
  chr <- as.character(seqnames(a))
  position <- c(start(a), end(a), start(b), end(b))
  dist <- abs(mean(c(start(a), end(a))) - start(promoters(b))[1])
  dist <- ceiling(dist/1000)
  left <- min(position) - 2000
  right <- max(position) + 2000
  return(paste(chr, " ", left, " ", right, " (", dist, "kb)", sep=""))
}

pdf(file="Correlated_eRNA_gene_pairs_GM.pdf", height=9, width=9)
par(mfrow=c(3,3))
my_palette <- colorRampPalette(c("blue", "red"), bias=1)(n = 9)
for(i in cluster_order){
  eRNA_list <- rownames(zscore)[k_cluster==i]
  for(j in 1:length(eRNA_list)){
    lab <- find_target(eRNA_list[j])
    if(length(lab)>0){
      for(k in 1:length(lab)){
        pvalue <- cor.test(e_mat[eRNA_list[j], ], g_mat[lab[k], ])$p.value
        region <- find_region(eRNA[mcols(eRNA)$id == eRNA_list[j]], gene[mcols(gene)$id %in% sym_2_ref(lab[k])])
        if(pvalue < 0.01){
          plot(e_mat[eRNA_list[j], ], g_mat[lab[k], ], col=my_palette, pch=16, cex=1.5,
               main=region, xlab=eRNA_list[j], ylab=lab[k])
          legend("topleft", paste("P <", signif(pvalue, 1)), bty="n")
        }
      }
    }
  }
}
dev.off()

###############################################################
# 5. Others
###############################################################
load("deg_matrix_GM_compare_with_30min.RData")
de_genes <- rownames(deg_matrix)[rowSums(deg_matrix[,1:8])>0]
tss <- promoters(gene, upstream = 1, downstream = 1)
flanking <- 500000
for(i in cluster_order){
  lab <- rownames(zscore)[k_cluster==i]
  tp <- eRNA[mcols(eRNA)$id %in% lab]
  start(tp) <- start(tp) - flanking
  end(tp) <- end(tp) + flanking
  hits <- findOverlaps(tp, tss)
  lab <- as.character(unique(mcols(tss[subjectHits(hits)])$id))
  lab <- ref_2_symbol(lab)
  lab <- lab[lab %in% de_genes]
  if(i==2){write.table(lab, quote=F, row.names = F, col.names = F, file="test.txt")}
  tp <- gene.rpkm.symbol[lab,]
  tp <- (tp - rowMeans(tp)) / apply(tp, 1, sd)
  print(length(lab))
  #boxplot(tp, ylab="Zscore")
  heatmap.2(tp, col=my_palette, Colv=F, Rowv=T, dendrogram="none",
            keysize=1.5, density.info="density",
            key.ylab="", key.xlab="", key.title="", trace="none")
}

par(mfrow=c(1,2))
my_palette <- colorRampPalette(c("blue", "red", "orange"), bias=1)(n = 9)
plot(eRNA.rpkm["MetaEnhancer_65_-", 2:10], gene.rpkm.symbol["IRF2",], col=my_palette, pch=16, cex=1.3,
     xlab="MetaEnhancer_65_-", ylab="IRF2")
plot(eRNA.rpkm["MetaEnhancer_250_-", 2:10], gene.rpkm.symbol["IRF1",], col=my_palette, pch=16, cex=1.3,
     xlab="MetaEnhancer_250_-", ylab="IRF1")

plot(eRNA.rpkm["MetaEnhancer_116_+", 2:10], gene.rpkm.symbol["TNFSF10",], col=my_palette, pch=16, cex=1.3,
     xlab="MetaEnhancer_116_+", ylab="TNFSF10")
