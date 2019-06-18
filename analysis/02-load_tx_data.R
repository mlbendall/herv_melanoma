#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

# Load samples
load(snakemake@input[["samp_rdata"]])

cat("--------- Sample table ---------\n")
samples %>% head

## Load transcript to gene mapping
ttg <- read.table(snakemake@input[["ttg_tsv"]], 
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F
                  )

cat("---------- Tx to Gene ----------\n")
ttg %>% head

## Load gene ID to symbol mapping
gsym <- read.table(snakemake@input[["gsym_tsv"]], 
                   sep = '\t',
                   header = T,
                   stringsAsFactors = F
                   )

cat("--------- Gene to Sym ----------\n")
gsym %>% head

## Load Kallisto counts
files <- snakemake@input[["h5_files"]]
names(files) <- sapply(files, function(x) strsplit(x, '/')[[1]][2])
files <- files[rownames(samples)]

txi <- tximport(files, type = "kallisto", tx2gene=ttg)

cat("--------- Count table ----------\n")
txi$counts[1:5, 1:5]

save(ttg, gsym, txi, file=snakemake@output[[1]])
