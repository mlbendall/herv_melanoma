#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)
library(DESeq2)

library(BiocParallel)
register(MulticoreParam(snakemake@threads))

#--- Parameters
padj_cutoff <- 0.1  # Adjusted p-value (FDR) cutoff
lfc_cutoff <- 1.0   # log2FoldChange (LFC) cutoff

# Contrasts
categories <- c('UVM.TP', 'SKCM.TP', 'SKCM.TM')
ctrst <- list(
    UM_v_CM.P = c('UVM.TP','SKCM.TP'), 
    UM_v_CM.M = c('UVM.TP','SKCM.TM'), 
    CM.P_v_CM.M = c('SKCM.TP','SKCM.TM')
)

# Use library size from alignment
use_libsize <- FALSE

#--- Load samples and metrics
load(snakemake@input[["samp_rdata"]])
load(snakemake@input[["metrics_rdata"]])

#--- Subset samples
sel <- samples$category %in% categories
sel.samples <- droplevels(samples[sel, ])
sel.metrics <- metrics[sel, ]

#--- Transcriptome DGE
load(snakemake@input[["tx_rdata"]])
sel.txi <- lapply(txi, function(x) if(is.matrix(x)) return(x[,sel]) else return(x))

stopifnot(all(colnames(txi$counts) == rownames(samples)))
stopifnot(all(colnames(sel.txi$counts) == rownames(sel.samples)))

dds.tx <- DESeqDataSetFromTximport(sel.txi, sel.samples, ~category)
if(use_libsize) sizeFactors(dds.tx) <- sel.metrics$ts.size_factor
dds.tx <- DESeq(dds.tx, parallel=T)

tform.tx <- varianceStabilizingTransformation(dds.tx, blind=FALSE)

results.tx <- lapply(ctrst, function(t)
    results(dds.tx, contrast=c("category", t[1], t[2]), alpha=padj_cutoff, parallel=T)
)
names(results.tx) <- names(ctrst)

cat("---------- tx result 1 ---------\n")
summary(results.tx[[1]])
cat("---------- tx result 2 ---------\n")
summary(results.tx[[2]])
cat("---------- tx result 3 ---------\n")
summary(results.tx[[3]])


#--- Load retrotranscriptome
load(snakemake@input[["rtx_rdata"]])
sel.rxi <- list(counts=rxi$counts[,sel])

stopifnot(all(colnames(rxi$counts) == rownames(samples)))
stopifnot(all(colnames(sel.rxi$counts) == rownames(sel.samples)))

dds.rtx <- DESeqDataSetFromMatrix(sel.rxi$counts, sel.samples, ~category)
if(use_libsize) sizeFactors(dds.rtx) <- sel.metrics$ts.size_factor
dds.rtx <- DESeq(dds.rtx, parallel=T)

tform.rtx <- varianceStabilizingTransformation(dds.rtx, blind=FALSE)

results.rtx <- lapply(ctrst, function(t)
    results(dds.rtx, contrast=c("category", t[1], t[2]), alpha=padj_cutoff, parallel=T)
)
names(results.rtx) <- names(ctrst)

cat("--------- rtx result 1 ---------\n")
summary(results.rtx[[1]])
cat("--------- rtx result 2 ---------\n")
summary(results.rtx[[2]])
cat("--------- rtx result 3 ---------\n")
summary(results.rtx[[3]])

save(ctrst, alpha, 
     dds.tx, tform.tx, results.tx,
     dds.rtx, tform.rtx, results.rtx,
     file=snakemake@output[[1]])
