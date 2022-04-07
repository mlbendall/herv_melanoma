#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

# load('analysis2/01-sample_data.Rdata')
load('analysis2/03-gene_data.Rdata')


samples.bap1 <- data.frame(
  sample_id =        c('SAMN13632402', 'SAMN13632403', 'SAMN13632404', 'SAMN13632405'),
  sample_sub =       c('A',            'B',            'C',            'D'),
  cell_line = factor(c('Mel202',       'Mel202',       'OMM2.3',       'OMM2.3')),
  bap1 =      factor(c('wt',           'KO',           'wt',           'KO'))
)

# Kallisto files
k_files <- file.path('add_samples', samples.bap1$sample_id, 'abundance.h5')
names(k_files) <- samples.bap1$sample_id
# Telescope files
t_files <- file.path('add_samples', samples.bap1$sample_id, 'telescope.report.tsv')
names(t_files) <- samples.bap1$sample_id
out_dat <- 'analysis2/31-bap1_counts.Rdata'


################################################################################
### Load Kallisto
################################################################################
txi <- tximport(k_files, type = "kallisto", tx2gene=ttg)

cat("------- Kallisto table ---------\n")
txi$counts[1:5, 1:4]

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
rxi$counts[1:5, 1:4]

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

# rename everything
counts.bap1.tx <- counts.tx
counts.bap1.rtx <- counts.rtx
counts.bap1.comb <- counts.comb
counts.bap1.herv <- counts.herv
counts.bap1.l1 <- counts.l1

save(counts.bap1.tx, counts.bap1.rtx, counts.bap1.comb, counts.bap1.herv, counts.bap1.l1, 
     file=out_dat)
rm(list=setdiff(ls(), c('counts.bap1.tx', 'counts.bap1.rtx', 'counts.bap1.comb', 'counts.bap1.herv', 'counts.l1')))
