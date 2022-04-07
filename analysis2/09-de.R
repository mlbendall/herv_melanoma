#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(UpSetR)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
  # load(snakemake@input[["deseq_rdata"]])
  # load(snakemake@input[["clin_rdata"]])
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load('analysis2/03-gene_data.Rdata')  
  load("analysis2/05-filter_counts.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  
  out_dat <- 'analysis2/09-de.Rdata'
  out_txt <- 'analysis2/09-de.txt'
}

p.cutoff <- 1e-3
lfc.cutoff <- 1.5

################################################################################
### HERV expression data
################################################################################
countDat <- mfilt.herv
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(samples)))
samples$clust.herv <- clust.df$clust.retro.k4

# Wald test
dds <- DESeq2::DESeqDataSetFromMatrix(countDat, samples, ~clust.herv + 0)
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

# Extract results
res <- list(
  "C1" = results(dds, contrast=c(+1, -1/3, -1/3, -1/3), alpha=p.cutoff),
  "C2" = results(dds, contrast=c(-1/3, +1, -1/3, -1/3), alpha=p.cutoff),
  "C3" = results(dds, contrast=c(-1/3, -1/3, +1, -1/3), alpha=p.cutoff),
  "C4" = results(dds, contrast=c(-1/3, -1/3, -1/3, +1), alpha=p.cutoff),
  "C12vC34" = results(dds, contrast=c(+1/2, +1/2, -1/2, -1/2), alpha=p.cutoff),
  "C1vC2"   = results(dds, contrast=c(  +1,   -1,    0,    0), alpha=p.cutoff),
  "C3vC4"   = results(dds, contrast=c(   0,    0,   +1,   -1), alpha=p.cutoff)
)

res <- lapply(res, function(r) {
  r$class <- loc_class[rownames(r),]$class
  r$display <- loc_class[rownames(r),]$display
  r
})

sig <- lapply(res, function(r) {
  s <- subset(r, padj < p.cutoff & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

sink(file=out_txt, split=T)
for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig[[n]])
}
sink()

save(dds, tform, res, sig, p.cutoff, lfc.cutoff, file=out_dat)
rm(list=setdiff(ls(), c('dds', 'tform', 'res', 'sig', 'p.cutoff', 'lfc.cutoff')))

