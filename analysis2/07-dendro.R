#! /usr/bin/env Rscript

library(tidyverse)
library(dendextend)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["clin_rdata"]])
  load(snakemake@input[["clust_rdata"]])
  # out_dat <- snakemake@output[["rdata"]]
  outdir <- snakemake@output[["outdir"]]
  
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  # out_dat <- 'analysis2/07-dendro.Rdata'
  outdir <- 'analysis2/07-dendro'
}

makePDF <- TRUE
if(makePDF) dir.create(outdir, showWarnings = F)
source('analysis2/constants.R')

# Use the tree for k=4
selK <- 4
cobj <- ccp.obj[[selK]]


if(exists("snakemake")) {
  # Rotate around the leftmost node
  xobj <- cobj$consensusTree %>% as.dendrogram()
  tmp <- xobj[[1]][[1]]
  xobj[[1]][[1]] <- xobj[[1]][[2]]
  xobj[[1]][[2]] <- tmp
  dobj <- xobj
  rm(tmp)
} else {
  # Interactive
  dobj <- click_rotate(cobj$consensusTree %>% as.dendrogram(), continue=T)
}

# Dendrogram ordering
dorder <- dobj %>% labels

mdata.o <- cbind(mdata, clust.df)
mdata.o <- mdata.o[dorder,]

if(makePDF) pdf(file.path(outdir, 'colorbars.pdf'), paper='letter')

par( oma = c(0,0,0,0), mgp = c(1,0.5,0), mar = c(24,4,2,2) )
dobj %>%
  dendextend::set('labels', rownames(mdata.o)) %>%
  dendextend::set('labels_cex', 0.4) %>% 
  dendextend::set('labels_col', cobj$clrs[[1]][dorder]) %>%
  plot(main='clust.retro.k4')

colbars <- data.frame(
  "EIF1AX" = my_pals$EIF1AX[as.character(mdata.o$EIF1AX)],
  "SRSF2" = my_pals$SRSF2[as.character(mdata.o$SRSF2)],
  "SF3B1" = my_pals$SF3B1[as.character(mdata.o$SF3B1)],
  "PLCB4" = my_pals$PLCB4[as.character(mdata.o$PLCB4)],
  "CYSLTR2" = my_pals$CYSLTR2[as.character(mdata.o$CYSLTR2)],
  "GNAQ" = my_pals$GNAQ[as.character(mdata.o$GNAQ)],
  "GNA11" = my_pals$GNA11[as.character(mdata.o$GNA11)],
  "BAP1" = my_pals$BAP1[mdata.o$BAP1],
  ".." = rep.int('#ffffff', nrow(mdata.o)),
  "chr8qISO" = my_pals$chr8qISO[as.character(mdata.o$chr8qISO)],
  "chr8qCN" = my_pals$chr8qCN[as.character(mdata.o$chr8qCN)],
  "chr3CN" = my_pals$chr3CN[as.character(mdata.o$chr3CN)],
  "metastasis" = my_pals$metastatic[mdata.o$metastatic],
  "death" = my_pals$death[mdata.o$death],
  "." = rep.int('#ffffff', nrow(mdata.o)),
  "SCNA" = my_pals$clust.SCNA[mdata.o$clust.SCNA],
  "DNA.methyl" = my_pals$clust.methyl[mdata.o$clust.methyl],
  "miRNA"    = my_pals$clust.miRNA[mdata.o$clust.miRNA],
  "PARADIGM" = my_pals$clust.paradigm[mdata.o$clust.paradigm],
  "lncRNA"   = my_pals$clust.lncRNA[mdata.o$clust.lncRNA],
  "mRNA"     = my_pals$clust.mRNA[mdata.o$clust.mRNA],
  "retro.k6" = my_pals$clust.retro.k6[mdata.o$clust.retro.k6],
  "retro.k5" = my_pals$clust.retro.k5[mdata.o$clust.retro.k5],
  "retro.k4" = my_pals$clust.retro.k4[mdata.o$clust.retro.k4]
)

  
dendextend::colored_bars(colbars, yshift=-.9)
if(makePDF) dev.off()



# 
# cstats <- fpc::cluster.stats(
#   d=dist(t(cDat), method = "euclidean"), 
#   clustering=as.numeric(clust.df$clust.retro.k4),
#   alt.clustering=as.numeric(mdata$clust.lncRNA)
# )
# cstats$corrected.rand
# cstats$vi
# 
# fisher.test(table(mdata$gender, mdata$clust.SCNA))
# fisher.test(table((ccp.obj[[4]]$consensusClass + 2) %% 4, ccp.obj[[4]]$consensusClass) )
# fisher.test(table((ccp.obj[[4]]$consensusClass + 2) %% 4, ccp.obj[[4]]$consensusClass) )
# 
# fisher.test(table(clust.df$clust.retro.k4, mdata$clust.lncRNA), simulate.p.value = TRUE)
# 
# 
# 
# fisher.test(table(mdata$PLCB4, mdata$clust.methyl))







