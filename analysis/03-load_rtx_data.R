#! /usr/bin/env Rscript

library(tidyverse)

#--- Load samples
load(snakemake@input[["samp_rdata"]])

cat("--------- Sample table ---------\n")
samples %>% head

#--- Load gene data
load(snakemake@input[["gene_rdata"]])

cat("----- Retroelement families ----\n")
retro.fam %>% head

cat("------- Retroelement table -----\n")
retro.annot %>% head

#--- Load Telescope counts
files <- snakemake@input[["tele_files"]]
names(files) <- sapply(files, function(x) strsplit(x, '/')[[1]][2])
files <- files[rownames(samples)]

cts <- lapply(1:length(files), 
              function(i) {
                  tmp <- read.table(files[i],
                                    sep='\t', header=T, stringsAsFactors=F)
                  ret <- data.frame(
                      transcript=retro.annot$locus, 
                      stringsAsFactors = F
                      ) %>%
                      left_join(tmp, by='transcript') %>%
                      select(gene_id=transcript, count=final_count)
                  
                  ret[is.na(ret)] <- 0
                  stopifnot(all(ret$gene_id == retro.annot$locus))
                  ret$gene_id <- NULL
                  names(ret) <- c(names(files)[i])
                  ret
              }) %>% bind_cols

row.names(cts) <- retro.annot$locus
rxi <- list(counts=cts)

cat("--------- Count table ----------\n")
rxi$counts[1:5, 1:5]

save(rxi, file=snakemake@output[[1]])
