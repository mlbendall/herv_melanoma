#! /usr/bin/env Rscript

library(tidyverse)
library(tximport)

# SF3B1 degron knockdown
samples.GSE66719<- read.table(file.path('add_samples', 'GSE66719.txt'), sep='\t')
samples.GSE66719$study <- 'GSE66719'

# Inhibition of Gaq and MEK
samples.GSE152705<- read.table(file.path('add_samples', 'GSE152705.txt'), sep='\t')
samples.GSE152705$study <- 'GSE152705'

# HDACi and mTOR
samples.GSE155452<- read.table(file.path('add_samples', 'GSE155452.txt'), sep='\t')
samples.GSE155452$study <- 'GSE155452'

# BAP1 wt vs. KO
samples.PRJNA596363 <- read.table(file.path('add_samples', 'PRJNA596363.txt'), sep='\t')
samples.PRJNA596363$study <- 'PRJNA596363'

# SF3B1 wt vs. mut
samples.GSE124720<- read.table(file.path('add_samples', 'GSE124720.txt'), sep='\t')
samples.GSE124720$study <- 'GSE124720'
samples.GSE124720 <- samples.GSE124720[samples.GSE124720$cancer=='UVM',]

# normal melanocyte study
samples.GSE176345<- read.table(file.path('add_samples', 'GSE176345.txt'), sep='\t')
samples.GSE176345$study <- 'GSE176345'

# BAP1 wt vs. mut
samples.GSE149920 <- read.table(file.path('add_samples', 'GSE149920.txt'), sep='\t')
samples.GSE149920$study <- 'GSE149920'

# UM patient sequencing 
samples.PRJNA386840<- read.table(file.path('add_samples', 'PRJNA386840.txt'), sep='\t')
samples.PRJNA386840$study <- 'PRJNA386840'

rbind(
  samples.GSE66719,    # Mel202
  samples.GSE152705,   # Mel202, OMM1.3
  samples.GSE155452,   # Mel202
  samples.PRJNA596363, # Mel202, OMM2.3
  samples.GSE124720,   # Mel270
  samples.GSE176345,   # 92.1, PIG1
  samples.GSE149920    # MM66, MP38, MP46, MP65, PDX4 
  #samples.PRJNA386840
) -> samples.add


load('analysis2/03-gene_data.Rdata')

# Kallisto files
k_files <- file.path('add_samples', samples.add$accession, 'abundance.h5')
names(k_files) <- samples.add$sample_id
# Telescope files
t_files <- file.path('add_samples', samples.add$accession, 'telescope.report.tsv')
names(t_files) <- samples.add$sample_id
out_dat <- 'analysis2/31-additional_counts.Rdata'


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
counts.add.tx <- counts.tx
counts.add.rtx <- counts.rtx
counts.add.comb <- counts.comb
counts.add.herv <- counts.herv
counts.add.l1 <- counts.l1

save(samples.add, counts.add.tx, counts.add.rtx, counts.add.comb, counts.add.herv, counts.add.l1, 
     file=out_dat)
rm(list=setdiff(ls(), c('samples.add', 'counts.add.tx', 'counts.add.rtx', 'counts.add.comb', 'counts.add.herv', 'counts.add.l1')))
