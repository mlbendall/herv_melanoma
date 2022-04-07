#! /usr/bin/env Rscript

library(tidyverse)
library(dendextend)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["clin_rdata"]])
  load(snakemake@input[["clust_rdata"]])
  out_pdf <- snakemake@output[["outpdf"]]
  out_dat <- snakemake@output[["rdata"]]
  
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  out_pdf <- 'analysis2/06-unsupervised/colorbars.pdf'
  out_dat <- 'analysis2/06a-dendro.Rdata'
}

makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)
source('analysis2/constants.R')

# Use the tree for k=4
selK <- 4
cobj <- ccp.obj[[selK]]

# Rotate around the leftmost node
if(exists("snakemake")) {
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

if(makePDF) pdf(out_pdf, paper='letter')

par( oma = c(0,0,0,0), mgp = c(1,0.5,0), mar = c(24,4,2,2) )
dobj %>%
  dendextend::set('labels', rownames(mdata.o)) %>%
  dendextend::set('labels_cex', 0.4) %>% 
  dendextend::set('labels_col', cobj$clrs[[1]][dorder]) %>%
  plot(main='HERV cluster')

colbars <- data.frame(
  "EIF1AX" = my_pals$EIF1AX[as.character(mdata.o$EIF1AX)],
  "SRSF2" = my_pals$SRSF2[as.character(mdata.o$SRSF2)],
  "SF3B1" = my_pals$SF3B1[as.character(mdata.o$SF3B1)],
  "PLCB4" = my_pals$PLCB4[as.character(mdata.o$PLCB4)],
  "CYSLTR2" = my_pals$CYSLTR2[as.character(mdata.o$CYSLTR2)],
  "GNAQ" = my_pals$GNAQ[as.character(mdata.o$GNAQ)],
  "GNA11" = my_pals$GNA11[as.character(mdata.o$GNA11)],
  "BAP1" = my_pals$BAP1[mdata.o$BAP1],
  "..-" = rep.int('#ffffff', nrow(mdata.o)),
  "chr8qISO" = my_pals$chr8qISO[as.character(mdata.o$chr8qISO)],
  "chr8qCN" = my_pals$chr8qCN[as.character(mdata.o$chr8qCN)],
  "chr3CN" = my_pals$chr3CN[as.character(mdata.o$chr3CN)],
  "D3/M3" = my_pals$D3M3[as.character(mdata.o$D3M3)],
  ".." = rep.int('#ffffff', nrow(mdata.o)),
  "stage" = my_pals$stage[as.character(mdata.o$stage)],
  "gender" = my_pals$gender[as.character(mdata.o$gender)],
  "metastasis" = my_pals$metastatic[mdata.o$metastatic],
  "death" = my_pals$death[mdata.o$death],
  "." = rep.int('#ffffff', nrow(mdata.o)),
  "SCNA" = my_pals$clust.SCNA[mdata.o$clust.SCNA],
  "DNA.methyl" = my_pals$clust.methyl[mdata.o$clust.methyl],
  "miRNA"    = my_pals$clust.miRNA[mdata.o$clust.miRNA],
  "PARADIGM" = my_pals$clust.paradigm[mdata.o$clust.paradigm],
  "lncRNA"   = my_pals$clust.lncRNA[mdata.o$clust.lncRNA],
  "mRNA"     = my_pals$clust.mRNA[mdata.o$clust.mRNA],
  # "retro.k6" = my_pals$clust.retro.k6[mdata.o$clust.retro.k6],
  # "retro.k5" = my_pals$clust.retro.k5[mdata.o$clust.retro.k5],
  ".x" = rep.int('#ffffff', nrow(mdata.o)),
  "HERV" = my_pals$clust.retro.k4[mdata.o$clust.retro.k4]
)

  
dendextend::colored_bars(colbars, yshift=-.9)

if(makePDF) dev.off()

save(dobj, dorder, mdata.o, file=out_dat)

################################################################################
### Immune landscape
################################################################################
celltypes <- c("B.cells.naive",
               "B.cells.memory",
               "T.cells.CD4.memory.activated",
               "T.cells.CD4.memory.resting",
               "T.cells.CD4.naive", 
               "T.cells.follicular.helper",
               "T.cells.CD8")
tmp  <- lapply(celltypes, function(n) {
  viridis::inferno(101)[floor(pmin(1, mdata.o[,n] / max(mdata.o[,n])) * 100) + 1]
})
names(tmp) <- celltypes

colbarsX <- cbind(
  data.frame(tmp),
  data.frame(
      "." = rep.int('#ffffff', nrow(mdata.o)),
      "leukfrac" = viridis::inferno(101)[floor(pmin(1, mdata.o$leukocyte_frac / 0.3) * 100) + 1],  
      ".x" = rep.int('#ffffff', nrow(mdata.o)),
      "HERV" = my_pals$clust.retro.k4[mdata.o$clust.retro.k4]
    )
)

if(makePDF) pdf('analysis2/06-unsupervised/colorbars_immune.pdf', paper='letter')

par( oma = c(0,0,0,0), mgp = c(1,0.5,0), mar = c(24,4,2,2) )
dobj %>%
  dendextend::set('labels', rownames(mdata.o)) %>%
  dendextend::set('labels_cex', 0.4) %>% 
  dendextend::set('labels_col', cobj$clrs[[1]][dorder]) %>%
  plot(main='HERV cluster')

dendextend::colored_bars(colbarsX, yshift=-.9)

if(makePDF) dev.off()