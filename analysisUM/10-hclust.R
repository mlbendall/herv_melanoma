#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(pvclust)
library(dendextend)

#--- Load DESeq objects
load(snakemake@input[["deseq_rdata"]])
load(snakemake@input[["clin_rdata"]])

load('analysisUM/05-deseq_int.Rdata')
load('analysisUM/01-load_clin_data.Rdata')


################################################################################
# Hierarchical clustering samples
################################################################################
set.seed(1234)

# Select most informative by median absolute deviation
tmat <- as.matrix(assay(tform.rtx))
mads <- apply(tmat, 1, mad)
tmat <- tmat[rev(order(mads))[1:2000],]

# Center by median
cmat <- sweep(tmat, 1, apply(tmat, 1, median, na.rm=T))

samples.pv <- pvclust::pvclust(cmat,
                               method.dist = 'correlation',
                               use.cor = "pairwise.complete.obs",
                               method.hclust = 'average', 
                               nboot=1000
                               )

pdf(snakemake@output[[1]], paper='letter')

colkey.clust <- c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6')

dorder <- samples.pv %>% as.dendrogram() %>% labels
mdata.dorder <- mdata[dorder,]
dlabels <- mdata.dorder$case_id
dcols <- colkey.clust[mdata.dorder$clust.SCNA]

par(mar=c(7,3,7,1))

samples.pv %>% 
    as.dendrogram() %>% 
    dendextend::set('labels', dlabels) %>%
    dendextend::set('labels_cex', 0.7) %>%
    dendextend::set('labels_col', dcols) %>%
    hang.dendrogram() %>%    
    plot(main="Cluster dendrogram with AU/BP values")

samples.pv %>% text(cex=0.7) 

dev.off()

