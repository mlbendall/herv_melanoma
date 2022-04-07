#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(ConsensusClusterPlus)
library(pheatmap)
# library(pvclust)
library(dendextend)

load('analysisUM/01-load_clin_data.Rdata')
load('analysisUM/05-deseq_int.Rdata')

# Select most informative by median absolute deviation
ntop <- 2000
tmat <- as.matrix(assay(tform.rtx))
mads <- apply(tmat, 1, mad)
tmat.top <- tmat[rev(order(mads))[1:ntop],]

# Center by median
cmat <- sweep(tmat.top, 1, apply(tmat.top, 1, median, na.rm=T))

# Spearman
title <- 'analysisUM/clust_results_S'
results.s <- ConsensusClusterPlus(cmat, maxK=20, reps=1000, pItem=0.8,
                                  pFeature=1, title=title, clusterAlg="hc",
                                  distance="spearman", seed=12345, plot="pdf")

icl.s <- calcICL(results.s, title=title, plot="pdf")

# Pearson
title <- 'analysisUM/clust_results_P'
results.p <- ConsensusClusterPlus(cmat, maxK=20, reps=1000, pItem=0.8,
                                  pFeature=1, title=title, clusterAlg="hc",
                                  distance="pearson", seed=12345, plot="pdf")

icl.p <- calcICL(results.p, title=title, plot="pdf")


# Based on delta area plot, the optimal number of clusters is 2
# We will use the spearman results
k <- 2
con.class <- results.s[[k]][['consensusClass']]
con.mat <- results.s[[k]][['consensusMatrix']]
con.tree <- results.s[[k]][['consensusTree']]

# Note: Two samples assigned to different clusters for spearman or pearson
# > con.class[con.class != results.p[[k]][['consensusClass']]]
# TCGA-VD-A8KK-01A TCGA-V4-A9EO-01A 
# (cluster 2 for spearman, 1 for pearson)
# > which(colnames(tmat) == 'TCGA-VD-A8KK-01A')
# > which(colnames(tmat) == 'TCGA-V4-A9EO-01A')

# Recalculate DESeq object
stopifnot(rownames(colData(dds.rtx)) == names(con.class))
colData(dds.rtx)$rtx_cluster <- factor(con.class, levels=c(1,2))
design(dds.rtx) <- ~rtx_cluster
dds.rtx <- DESeq(dds.rtx, parallel=T)

# Criteria
alpha <- 0.01
lfc <- 1
ntop.hm <- 50

res <- results(dds.rtx, contrast=c("rtx_cluster", 1 ,2), alpha=alpha)
sig <- subset(res, padj<alpha & abs(log2FoldChange)>=lfc)

nameorder <- rownames(sig[order(sig$padj),])

# gets only HERV
topgenes <- nameorder[grep('^L1', nameorder,invert=T)][1:ntop.hm]
mat <- assay(tform.rtx)[topgenes,]

stopifnot(all(rownames(mdata) == colnames(mat)))
rowdist <- as.dist((1 - cor(t(mat), method='spearman'))/2)
pheatmap(mat, 
         scale="row", 
         clustering_distance_rows = rowdist,
         clustering_method="average",
         annotation_col=mdata[,c('clust.SCNA', 'clust.methyl',  'clust.mRNA', 'metastatic', 'BAP1')]
)



1 - cor(t(mat), method='spearman')


