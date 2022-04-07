#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(UpSetR)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
} else {
  cat("Not using snakemake\n")
  load('analysis2/09-de.Rdata')
  out_pdf <- 'analysis2/09-de/upset.pdf'
  out_dat <- 'analysis2/09a-upset.Rdata'
}


makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)

# Activated HERVs
upvars <- lapply(sig[1:4], function(r) rownames(subset(r, log2FoldChange>0)))
# Silenced HERVs
downvars <- lapply(sig[1:4], function(r) rownames(subset(r, log2FoldChange<0)))


if(makePDF) pdf(out_pdf, wi onefile = T)
upset(fromList(upvars), sets=c("C4", "C3", "C2", "C1"),  
      keep.order = T, order.by='degree', decreasing=F)

upset(fromList(downvars), sets=c("C4", "C3", "C2", "C1"),  
      keep.order = T, order.by='degree', decreasing=F)

# intersections=list(list("C1"),list("C2"),list("C3"),list("C4"), 
#                    list("C1","C2"), list("C1","C3"),list("C1","C4"),list("C2","C3"),list("C2","C4"), list("C3","C4")
# ),
if(makePDF) dev.off()

up.binmat <- fromList(upvars)
rn <- do.call(c, upvars)
rn <- rn[!duplicated(rn)]
rownames(up.binmat) <- rn
rm(rn)

dn.binmat <- fromList(downvars)
rn <- do.call(c, downvars)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat) <- rn
rm(rn)



save(up.binmat, dn.binmat, file=out_dat)
rm(list=setdiff(ls(), c('up.binmat', 'dn.binmat')))
