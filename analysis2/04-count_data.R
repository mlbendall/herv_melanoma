#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["gene_rdata"]])
  # Kallisto files
  k_files <- snakemake@input[["h5_files"]]
  names(k_files) <- sapply(k_files, function(x) strsplit(x, '/')[[1]][2])
  k_files <- k_files[rownames(samples)]
  # Telescope files
  t_files <- snakemake@input[["tele_files"]]
  names(t_files) <- sapply(t_files, function(x) strsplit(x, '/')[[1]][2])
  t_files <- t_files[rownames(samples)]
  out_dat <- snakemake@output[[1]]
} else {
  cat("Not using snakemake\n")
  load('analysis2/01-sample_data.Rdata')
  load('analysis2/03-gene_data.Rdata')
  # Kallisto files
  k_files <- file.path('samples', samples$sample_id, 'abundance.h5')
  names(k_files) <- samples$sample_id
  # Telescope files
  t_files <- file.path('samples', samples$sample_id, 'telescope.report.tsv')
  names(t_files) <- samples$sample_id
  out_dat <- 'analysis2/04-count_data.Rdata'
}

cat("--------- Sample table ---------\n")
samples %>% head

cat("---------- Tx to Gene ----------\n")
ttg %>% head

cat("--------- Gene to Sym ----------\n")
gsym %>% head

cat("----- Retroelement families ----\n")
retro.fam %>% head

cat("------- Retroelement table -----\n")
retro.annot %>% head

################################################################################
### Load Kallisto
################################################################################
txi <- tximport(k_files, type = "kallisto", tx2gene=ttg)

cat("------- Kallisto table ---------\n")
txi$counts[1:5, 1:5]


################################################################################
### Load Telescope
################################################################################
cts <- lapply(1:length(t_files), 
              function(i) {
                tmp <- read.table(t_files[i],
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
                names(ret) <- c(names(t_files)[i])
                ret
              }) %>% bind_cols

row.names(cts) <- retro.annot$locus
rxi <- list(counts=cts)

cat("------- Telescope table --------\n")
rxi$counts[1:5, 1:5]

################################################################################
### Combined counts
################################################################################
stopifnot(all(colnames(txi$counts) == colnames(rxi$counts)))

counts.tx <- round(txi$counts)
counts.rtx <- rxi$counts
counts.comb <- rbind(counts.tx, counts.rtx)

stopifnot(all(rownames(counts.rtx) == rownames(retro.annot)))
counts.herv <- counts.rtx[retro.annot$class == 'HERV',]
counts.l1 <- counts.rtx[retro.annot$class == 'L1',]

cat(sprintf("%d genes\n", nrow(counts.tx)))
cat(sprintf("%d retro genes\n", nrow(counts.rtx)))
cat(sprintf("%d combined genes\n", nrow(counts.comb)))
cat(sprintf("%d HERV\n", nrow(counts.herv)))
cat(sprintf("%d L1\n", nrow(counts.l1)))

save(counts.tx, counts.rtx, counts.comb, counts.herv, counts.l1, file=out_dat)
rm(list=setdiff(ls(), c('counts.tx', 'counts.rtx', 'counts.comb', 'counts.herv', 'counts.l1')))
