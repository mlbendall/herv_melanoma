#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["count_rdata"]])
  out_dat <- snakemake@output[[1]]
} else {
  cat("Not using snakemake\n")
  load('analysis2/04-count_data.Rdata')
  out_dat <- 'analysis2/05-filter_counts.Rdata'
}

cutoff.count <- 5
cutoff.samp <- floor(ncol(counts.comb) * 0.05)

################################################################################
### Minimal filter: by total count
################################################################################
cat(sprintf("Minimal filtering: Keep genes with more than %d total counts.\n", cutoff.count))

mfilt.tx <- counts.tx[rowSums(counts.tx) > cutoff.count, ]
mfilt.rtx <- counts.rtx[rowSums(counts.rtx) > cutoff.count, ]
mfilt.comb <- rbind(mfilt.tx, mfilt.rtx)
mfilt.herv <- counts.herv[rowSums(counts.herv) > cutoff.count, ]
mfilt.l1 <- counts.herv[rowSums(counts.l1) > cutoff.count, ]

cat(sprintf("%d genes retained (%d original)\n", nrow(mfilt.tx), nrow(counts.tx)))
cat(sprintf("%d retro genes retained (%d original)\n", nrow(mfilt.rtx), nrow(counts.rtx)))
cat(sprintf("%d combined genes retained (%d original)\n", nrow(mfilt.comb), nrow(counts.comb)))
cat(sprintf("%d HERVs retained (%d original)\n", nrow(mfilt.herv), nrow(counts.herv)))
cat(sprintf("%d L1s retained (%d original)\n", nrow(mfilt.l1), nrow(counts.l1)))


################################################################################
### Moderate filter: by count and sample
################################################################################
cat(sprintf("Filtering: Keep genes with more than %d counts in more than %d samples.\n", cutoff.count, cutoff.samp))

filt.tx <- counts.tx[rowSums(counts.tx > cutoff.count) > cutoff.samp, ]
filt.rtx <- counts.rtx[rowSums(counts.rtx > cutoff.count) > cutoff.samp, ]
filt.comb <- rbind(filt.tx, filt.rtx)
filt.herv <- counts.herv[rowSums(counts.herv > cutoff.count) > cutoff.samp, ]
filt.l1 <- counts.herv[rowSums(counts.l1 > cutoff.count) >cutoff.samp, ]

cat(sprintf("%d genes retained (%d original)\n", nrow(filt.tx), nrow(counts.tx)))
cat(sprintf("%d retro genes retained (%d original)\n", nrow(filt.rtx), nrow(counts.rtx)))
cat(sprintf("%d combined genes retained (%d original)\n", nrow(filt.comb), nrow(counts.comb)))
cat(sprintf("%d HERVs retained (%d original)\n", nrow(filt.herv), nrow(counts.herv)))
cat(sprintf("%d L1s retained (%d original)\n", nrow(filt.l1), nrow(counts.l1)))

save(mfilt.tx, mfilt.rtx, mfilt.comb, mfilt.herv, mfilt.l1, 
     filt.tx,  filt.rtx,  filt.comb,  filt.herv,  filt.l1, 
     file=out_dat)

rm(list=setdiff(ls(), c('mfilt.tx', 'mfilt.rtx', 'mfilt.comb', 'mfilt.herv', 'mfilt.l1',
                        'filt.tx',  'filt.rtx',  'filt.comb',  'filt.herv',  'filt.l1')))
