#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

#--- Load samples
load(snakemake@input[["samp_rdata"]])

cat("--------- Sample table ---------\n")
samples %>% head

#--- Load gene data
load(snakemake@input[["gene_rdata"]])

cat("---------- Tx to Gene ----------\n")
ttg %>% head

cat("--------- Gene to Sym ----------\n")
gsym %>% head

#--- Load Kallisto counts
files <- snakemake@input[["h5_files"]]
names(files) <- sapply(files, function(x) strsplit(x, '/')[[1]][2])
files <- files[rownames(samples)]

txi <- tximport(files, type = "kallisto", tx2gene=ttg)

cat("--------- Count table ----------\n")
txi$counts[1:5, 1:5]

save(txi, file=snakemake@output[[1]])
