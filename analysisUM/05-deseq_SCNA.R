#! /usr/bin/env R

library(tidyverse)
library(tximport)
library(DESeq2)

library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#--- Parameters
padj_cutoff <- 0.1  # Adjusted p-value (FDR) cutoff
lfc_cutoff <- 1.0   # log2FoldChange (LFC) cutoff

# Contrasts
categories <- c('UVM.TP')
ctrst <- list(
    C1_v_C2 = c('C2','C1'),
    C1_v_C3 = c('C3','C1'),
    C1_v_C4 = c('C4','C1'),
    C2_v_C3 = c('C3','C2'),
    C2_v_C4 = c('C4','C2'),
    C3_v_C4 = c('C4','C3') 
)

# Use library size from alignment
use_libsize <- FALSE

#--- Load samples and metrics
load(snakemake@input[["clin_rdata"]])
load(snakemake@input[["metrics_rdata"]])

#--- Subset samples
sel <- mdata$category %in% categories
sel.mdata <- droplevels(mdata[sel, ])
sel.metrics <- metrics[sel, ]

#--- Transcriptome DGE
load(snakemake@input[["tx_rdata"]])
sel.txi <- lapply(txi, function(x) if(is.matrix(x)) return(x[,sel]) else return(x))

stopifnot(all(colnames(txi$counts) == rownames(mdata)))
stopifnot(all(colnames(sel.txi$counts) == rownames(sel.mdata)))

dds.tx <- DESeqDataSetFromTximport(sel.txi, sel.mdata, ~clust.SCNA)
if(use_libsize) sizeFactors(dds.tx) <- sel.metrics$ts.size_factor
dds.tx <- DESeq(dds.tx, parallel=T)

tform.tx <- varianceStabilizingTransformation(dds.tx, blind=FALSE)

#--- Load retrotranscriptome
load(snakemake@input[["rtx_rdata"]])
sel.rxi <- list(counts=rxi$counts[,sel])

stopifnot(all(colnames(rxi$counts) == rownames(mdata)))
stopifnot(all(colnames(sel.rxi$counts) == rownames(sel.mdata)))

dds.rtx <- DESeqDataSetFromMatrix(sel.rxi$counts, sel.mdata, ~clust.SCNA)
if(use_libsize) sizeFactors(dds.rtx) <- sel.metrics$ts.size_factor
dds.rtx <- DESeq(dds.rtx, parallel=T)

tform.rtx <- varianceStabilizingTransformation(dds.rtx, blind=FALSE)



results.rtx <- lapply(ctrst, function(t)
    results(dds.rtx, contrast=c("clust.SCNA", t[1], t[2]), alpha=padj_cutoff, parallel=T)
)



save(dds.tx, tform.tx,
     dds.rtx, tform.rtx,
     file=snakemake@output[[1]])
