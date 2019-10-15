setwd("~/Documents/IFN_enhancer/R/")
library(RColorBrewer)
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(rtracklayer)
library(cluster)
library(vioplot)
source("Script/ref_2_symbol_matrix.R")
source("Script/sym_2_ref.R")
source("Script/liftover_function.R")
load("Data/gene_expression_GM_unfiltered.RData")
colnames(gene.rpkm) <- substr(colnames(gene.rpkm), 3, 1000L)

load("~/Documents/IFN_enhancer/R/Data/ref2sym_dedup.RData")
ref_2_symbol_local <- function(x){
  # stopifnot(all(x %in% ref2sym[,1]))
  # tp <- 1:nrow(ref2sym)
  # names(tp) <- ref2sym[,1]
  # return(ref2sym[tp[x], 2])
  output <- NULL
  for(i in 1:length(x)){
    if(x[i] %in% ref2sym[,1]){output <- c(output, ref2sym[which(ref2sym[,1]==x[i])[1], 2])}
    else{output <- c(output, NA)}
  }
  return(output)
}

fold_change <- function(x, y, pseudo=0.1){
  log_fold <- log2((x+pseudo)/(y+pseudo))
  return(log_fold)
}
amplitude_index <- function(x, abs=T){
  fc <- log2((x[,-1]+0.1) /(x[,1]+0.1))
  if(abs){return(apply(abs(fc), 1, max))}
  else{
    lab <- apply(abs(fc), 1, function(x){xx<-which(x==max(x));return(xx[1])})
    tp <- NULL
    for(i in 1:length(lab)){
      tp <- c(tp, fc[i, lab[i]])
    }
    names(tp) <- rownames(x)
    return(tp)
  }
}
continuity_index <- function(x){
  return(diag(cor(t(x[, -ncol(x)]), t(x[, -1]), method="spearman")))
}
find_pair <- function(a, b, a_mat, b_mat, distance, select_best=T, output_symbol=T){
  # Make sure that records of a/b matches with rows of a/b_mat
  id_a <- NULL
  id_b <- NULL
  id_a_all <- rownames(a_mat)
  id_b_all <- rownames(b_mat)
  a <- promoters(a, upstream=1, downstream=0)
  b <- promoters(b, upstream=1, downstream=0)
  for(i in 1:length(a)){
    inner <- flank(a[i], width=distance[1], both=T)
    outer <- flank(a[i], width=distance[2], both=T)  # Region flanking TSS
    tp <- setdiff(outer, inner)
    hits <- unique(subjectHits(findOverlaps(tp, b, ignore.strand=T)))
    if(length(hits) == 0){next}
    if(select_best){
      lab <- order(cor(t(b_mat)[,hits], a_mat[i,]), decreasing = T)[1]  # Incase of multiple enhancers, select the one that is the most similar to the gene
      hits <- hits[lab]
    }
    tp <- id_a_all[i]
    if(output_symbol){tp <- ref_2_symbol(tp)}
    if(length(tp) != 1){next}
    tp <- rep(tp, length(hits))
    id_a <- c(id_a, tp)
    tp <- id_b_all[hits]
    id_b <- c(id_b, tp)
  }
  return(data.frame(id_a, id_b, stringsAsFactors = F))
}
GRange_index <- function(region, id){
  # call records from GRange objects by the order of ids provided
  tp <- c(1:length(region))
  names(tp) <- as.character(mcols(region)$id)
  return(region[tp[id]])
}
PE_plot <- function(gid, eid, distance, col, black_list=NULL){
  par(mfrow=c(2,1), mar=c(3,4,1,1))
  rg <- GRange_index(gene, gid)
  mg <- gene.rpkm[gid, -1]
  re <- GRange_index(eRNA, eid)
  me <- eRNA.rpkm[eid, -1]
  PE <- find_pair(a=rg, b=re, a_mat=mg, b_mat=me, distance=distance)

  # Plot of gene
  lab <- unique(PE[,1])
  lab <- lab[!(lab %in% black_list)]
  tp <- fold_change(mg[lab,-1], mg[lab,1])
  plot(-1,-1, xlim=c(0, 11), ylim=c(-3,4), xaxt="n", ylab="Gene Log2 Fold-change", xlab="")
  axis(side=1, at=c(1:10), labels=colnames(tp))
  for(j in c(1:10)){vioplot(tp[,j], at=j, add=T, col=col[1])}
  text(1, 2.8, paste("n=", length(lab), sep=""))
  abline(h=0, lty="dashed", col="grey", lwd=1.5)

  # Plot of eRNA
  lab <- unique(PE[,2])
  lab <- lab[!(lab %in% black_list)]
  tp <- fold_change(me[lab,-1], me[lab,1])
  plot(-1,-1, xlim=c(0, 11), ylim=c(-3,4), xaxt="n", ylab="eRNA Log2 Fold-change")
  axis(side=1, at=c(1:10), labels=colnames(tp))
  for(j in c(1:10)){vioplot(tp[,j], at=j, add=T, col=col[2])}
  text(1, 2.8, paste("n=", length(lab), sep=""))
  abline(h=0, lty="dashed", col="grey", lwd=1.5)
  return(lab)
}
pair_region <- function(a, b){
  a <- GRange_index(gene, a)
  b <- GRange_index(eRNA, b)
  tp <- c(start(a), end(a), start(b), end(b))
  tp <- paste(seqnames(a), min(tp), max(tp))
  return(tp)
}
################################################################
# 1. Load data
################################################################

# 1.1 DE genes
# Criteria:
# 1. max(RPKM) > 1
# 2. Amplitude index > 1
# 3. Continuity index > 0.2
load(file="Data/gene_expression_GM_unfiltered.RData")
colnames(gene.rpkm) <- substr(colnames(gene.rpkm), 3, 1000L)
rpkm <- gene.rpkm[,-1]
rpkm <- rpkm[apply(rpkm, 1, max)>1, ] # 8484 genes remained
ai_rpkm <- amplitude_index(rpkm[,1:9], abs=F)   # Only consider early time points
ci_rpkm <- continuity_index(rpkm)
index_rpkm <- ai_rpkm * ci_rpkm

deg <- names(index_rpkm)[(ci_rpkm>0.2 & ai_rpkm>log2(2))]
deg.data <- cbind(rpkm, ai_rpkm, ci_rpkm, index_rpkm)[deg, ]
tp <- c(1:length(gene))
names(tp) <- as.character(gene$id)
deg.region <- gene[tp[deg]]  # Select corresponding regions

# 1.2 DE enhancers
# Criteria:
# 1. max(RPKM) > 0.5
# 2. Amplitude index > 1
# 3. Continuity index > 0.2
load(file="Data/eRNA_expression_GM_unfiltered.RData")
colnames(eRNA.rpkm) <- substr(colnames(eRNA.rpkm), 3, 1000L)
eRNA_rpkm <- eRNA.rpkm[,-1]
eRNA_rpkm <- eRNA_rpkm[apply(eRNA_rpkm, 1, max)>0.5, ]
ai_eRNA <- amplitude_index(eRNA_rpkm[,1:9], abs=F)
ci_eRNA <- continuity_index(eRNA_rpkm)
index_eRNA <- ai_eRNA*ci_eRNA

dee <- names(index_eRNA)[ci_eRNA>0.2 & ai_eRNA>log2(2)]
dee.data <- cbind(eRNA_rpkm, ai_eRNA, ci_eRNA, index_eRNA)[dee, ]
tp <- c(1:length(eRNA))
names(tp) <- as.character(eRNA$id)
enh.region <- eRNA[tp[dee]]

#ind_eRNA <- names(ci_eRNA)[ci_eRNA>0.2 & amplitude_index(eRNA_rpkm, abs=F)>log2(2)]

################################################################
# 2. DE gene classification
################################################################

heat.deg <- log2((gene.rpkm[deg, -c(1,2)]+0.1) / (gene.rpkm[deg,2]+0.1))
# Data clustering
dis_matrix <- (1-cor(t(heat.deg), method="pearson"))/2
set.seed(1)
centers <- 2
cluster_order <- c(1,2)
kmedoids <- pam(dis_matrix, k=centers, diss=T)  # Method: Partitioning Around Medoids
k_cluster <- kmedoids$clustering

# Reordering
ind <- cluster_order[k_cluster] * 10^8
for(i in 1:centers){
  lab <- which(k_cluster==i)
  ind[lab] <- ind[lab] - apply(heat.deg[lab, ], 1, max)
}
od <- order(ind, decreasing = F)
heat.deg <- heat.deg[od, ]
k_cluster <- k_cluster[od]

# Plotting
pdf("Figure/Diff_gene_heatmap.pdf")
my_palette <- my_palette <- colorRampPalette(c("blue", "white", "red"), bias=1)(n = 20)
heatmap.2(heat.deg, col=my_palette, Colv=F, Rowv=F, dendrogram="none",
          keysize=1.5, density.info="density", labRow=F,
          key.ylab="", key.xlab="", key.title="", trace="none")
dev.off()

write.table(ref_2_symbol(deg), "Data/Inducible_genes.symbol.list", quote=F, row.names=F, col.names=F)
write.table(deg, "Data/Inducible_genes.list", quote=F, row.names=F, col.names=F)

################################################################
# 3. DE eRNA classification
################################################################

heat.dee <- log2((eRNA.rpkm[dee, -c(1,2)]+0.1) / (eRNA.rpkm[dee,2]+0.1))
# Data clustering
dis_matrix <- (1-cor(t(heat.dee), method="pearson"))/2
set.seed(1)
centers <- 2
cluster_order <- c(2,1) # Cluster cluster_order[1]/[2]/[3] as: Top/mid/bot
kmedoids <- pam(dis_matrix, k=centers, diss=T)  # Method: Partitioning Around Medoids
k_cluster <- kmedoids$clustering

# Reordering
ind <- cluster_order[k_cluster] * 10^8
for(i in 1:centers){
  lab <- which(k_cluster==i)
  ind[lab] <- ind[lab] - apply(heat.dee[lab, ], 1, max)
}
od <- order(ind, decreasing = F)
heat.dee <- heat.dee[od, ]
k_cluster <- k_cluster[od]

# Plotting
pdf("Figure/Diff_eRNA_heatmap.pdf")
my_palette <- my_palette <- colorRampPalette(c("blue", "white", "red"), bias=1)(n = 20)
heatmap.2(heat.dee, col=my_palette, Colv=F, Rowv=F, dendrogram="none",
          keysize=1.5, density.info="density", labRow=F,
          key.ylab="", key.xlab="", key.title="", trace="none")
dev.off()

write.table(dee, "Data/Inducible_eRNA.list", quote=F, row.names=F, col.names=F)

################################################################
# 4. Associate inducible enhancers to inducible genes
################################################################

# Now we decided to use 200kb as the cutoff
selected_eRNA <- dee
selected_gene <- deg
ind_e <- GRange_index(eRNA, selected_eRNA)
ind_g <- GRange_index(gene, selected_gene)
mat_e <- dee.data[selected_eRNA, c(1:11)]
mat_g <- deg.data[selected_gene, c(1:11)]
PE <- find_pair(a=ind_g, b=ind_e, a_mat=mat_g, b_mat=mat_e, distance=c(0, 2*10^5), select_best = F)

write.table(PE, "Data/PE.list", quote=F, row.names=F, col.names=F, sep="\t")

# Put the inducible EP pairs and related information into a table

tp <- find_pair(a=ind_g, b=ind_e, a_mat=mat_g, b_mat=mat_e, distance=c(0, 2*10^5), select_best = F, output_symbol = F)
tp$symbol <- ref_2_symbol_local(tp[,1])
tp$coord <- "chr0"
tp$enh_hg18 <- "chr0"
tp$enh_hg19 <- "chr0"
for(i in 1:nrow(tp)){
  tp$coord[i] <- pair_region(tp[i,1], tp[i,2])
  tp_hg18 <- GRange_index(eRNA, tp[i,2])
  tp_hg19 <- liftOver_function(tp_hg18)[1]
  tp$enh_hg18[i] <- paste(seqnames(tp_hg18), start(tp_hg18), end(tp_hg18))
  tp$enh_hg19[i] <- paste(seqnames(tp_hg19), start(tp_hg19), end(tp_hg19))
}

a <- promoters(GRange_index(gene, tp[,1]), 0,0)
b <- promoters(GRange_index(eRNA, tp[,2]), 0,0)
tp$EP_d <- round((start(b)-start(a))/1000)
lab <- which(strand(a)=="-")
tp$EP_d[lab] <- tp$EP_d[lab] * (-1)
tp$EP_d <- paste(tp$EP_d, "kb", sep="")
tp <- tp[order(tp$symbol),]
tp <- cbind(tp, rpkm[tp$id_a, ])
tp <- cbind(tp, eRNA_rpkm[tp$id_b, ])
write.csv(tp[,], "Data/Induced_EP.csv", row.names = F)

#################################################################################################

# For enhancers without inducible target within 200kb, increase search space to 1Mb
selected_eRNA <- dee[-which(dee %in% PE$id_b)]
selected_gene <- deg
ind_e <- GRange_index(eRNA, selected_eRNA)
ind_g <- GRange_index(gene, selected_gene)
mat_e <- dee.data[selected_eRNA, c(1:11)]
mat_g <- deg.data[selected_gene, c(1:11)]
PE_others <- find_pair(a=ind_g, b=ind_e, a_mat=mat_g, b_mat=mat_e, distance=c(2*10^5, 1*10^6), select_best = F)

write.table(PE_others, "Data/PE_1MB.list", quote=F, row.names=F, col.names=F, sep="\t")

# Put the inducible EP pairs and related information into a table

tp <- find_pair(a=ind_g, b=ind_e, a_mat=mat_g, b_mat=mat_e, distance=c(2*10^5, 1*10^6), select_best = F, output_symbol = F)
tp$symbol <- ref_2_symbol_local(tp[,1])
tp$coord <- "chr0"
tp$enh_hg18 <- "chr0"
tp$enh_hg19 <- "chr0"
for(i in 1:nrow(tp)){
  tp$coord[i] <- pair_region(tp[i,1], tp[i,2])
  tp_hg18 <- GRange_index(eRNA, tp[i,2])
  tp_hg19 <- liftOver_function(tp_hg18)[1]
  tp$enh_hg18[i] <- paste(seqnames(tp_hg18), start(tp_hg18), end(tp_hg18))
  tp$enh_hg19[i] <- paste(seqnames(tp_hg19), start(tp_hg19), end(tp_hg19))
}

a <- promoters(GRange_index(gene, tp[,1]), 0,0)
b <- promoters(GRange_index(eRNA, tp[,2]), 0,0)
tp$EP_d <- round((start(b)-start(a))/1000)
lab <- which(strand(a)=="-")
tp$EP_d[lab] <- tp$EP_d[lab] * (-1)
tp$EP_d <- paste(tp$EP_d, "kb", sep="")
tp <- tp[order(tp$symbol),]
tp <- cbind(tp, rpkm[tp$id_a, ])
write.csv(tp[,], "Data/Induced_EP_1MB.csv", row.names = F)

################################################################
# 5. Number of inducible genes decreases with distance to enhancers
################################################################

step = 1*10^5
nbin = 10

# 5.1 Foreground
freq = NULL
selected_eRNA <- dee
selected_gene <- deg
ind_e <- GRange_index(eRNA, selected_eRNA)
ind_g <- GRange_index(gene, selected_gene)
mat_e <- dee.data[selected_eRNA, c(1:11)]
mat_g <- deg.data[selected_gene, c(1:11)]

for(i in 1:nbin){
  tp <- find_pair(a=ind_g, b=ind_e, a_mat=mat_g, b_mat=mat_e, distance=c((i-1)*step, i*step), select_best = F, output_symbol = F)
  freq = c(freq, nrow(tp))
}

# 5.2 Background
freq_bg = NULL
selected_eRNA <- dee
selected_gene <- rownames(gene.rpkm)
ind_e <- GRange_index(eRNA, selected_eRNA)
ind_g <- GRange_index(gene, selected_gene)
mat_e <- dee.data[selected_eRNA, c(1:11)]
mat_g <- gene.rpkm[selected_gene, c(1:11)]

for(i in 1:nbin){
  tp <- find_pair(a=ind_e, b=ind_g, a_mat=mat_e, b_mat=mat_g, distance=c((i-1)*step, i*step), select_best = F, output_symbol = F)
  freq_bg = c(freq_bg, nrow(tp))
}

# 5.3 Plot
fg = freq / max(freq)
bg = freq_bg / max(freq_bg)

pdf("Figure/Inducible_genes_near_enhancer.pdf")
postscript("Figure/Inducible_genes_near_enhancer.eps")
plot((c(1:nbin)-0.5)*step/10^6, fg, xlab="Distance to Enhancers (MB)", ylab="Normalized Target Occurence", main="",
     xaxp=c(0, nbin*step/10^6, nbin), yaxp=c(0, 1.1, 11), xlim=c(0, nbin*step/10^6), ylim=c(0, 1.2),
     col="blue", cex=1.5)
points((c(1:nbin)-0.5)*step/10^6, bg, col="grey", cex=1.5)
grid(ny=NULL, nx=NA, lwd=1.5)
legend("top", c("All Genes", "Inducible Genes"), ncol=2, col=c("grey", "blue"), cex=1)
dev.off()
