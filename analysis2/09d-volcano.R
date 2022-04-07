#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(gridExtra)


#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
  # load(snakemake@input[["deseq_rdata"]])
  # load(snakemake@input[["clin_rdata"]])
} else {
  cat("Not using snakemake\n")
  # load("analysis2/01-sample_data.Rdata")  
  # load("analysis2/02-clinical_data.Rdata")
  # load('analysis2/03-gene_data.Rdata')  
  # load("analysis2/05-filter_counts.Rdata")
  # load("analysis2/06-unsupervised.Rdata")
  load('analysis2/03-gene_data.Rdata')
  load('analysis2/06a-dendro.Rdata')
  load('analysis2/09-de.Rdata')
  load('analysis2/09a-upset.Rdata')
  out_pdf1 <- 'analysis2/09-de/volcano1.pdf'
  out_pdf2 <- 'analysis2/09-de/volcano2.pdf'
}


makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)

source('analysis2/constants.R')

pl <- lapply(names(res)[1:4], function(n) {
  r <- res[[n]]
  rs <- sig[[n]]
  rs <- subset(rs, class=='HERV')
  
  sellab <- unique(c(
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% filter(log2FoldChange>0) %>% slice_min(padj, n=3))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% filter(log2FoldChange<0) %>% slice_min(padj, n=3))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% slice_min(log2FoldChange, n=3))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% slice_max(log2FoldChange, n=3))$loc
  ))
  pointsize <- ifelse(rownames(r) %in% sellab, 2.5, 1.0)
  keyvals.shape <- ifelse(r$class=='HERV', 16, 1)
  names(keyvals.shape)[keyvals.shape==16] <- 'HERV'
  names(keyvals.shape)[keyvals.shape==1] <- 'other'
  
  EnhancedVolcano(r, 
                  lab=r$display, selectLab = sellab,
                  xlim=c(-6,6),
                  ylim=c(0,40),
                  pointSize = pointsize,
                  drawConnectors = T,
                  widthConnectors = 0.5,
                  col = c("grey30", "grey30", "grey30", "red2"),
                  x = 'log2FoldChange', y = 'padj',
                  labSize = 3.5,
                  axisLabSize = 8,
                  shapeCustom = keyvals.shape,
                  title = n, subtitle=NULL, 
                  legendPosition='none',
                  xlab=NULL, ylab=NULL, caption = NULL,
                  pCutoff=p.cutoff, FCcutoff=lfc.cutoff) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
})

if(makePDF) pdf(out_pdf1, width=16, height=4)
grid.arrange(grobs=pl, nrow=1)
if(makePDF) dev.off()


pl2 <- lapply(names(res), function(n) {
  r <- res[[n]]
  rs <- sig[[n]]
  rs <- subset(rs, class=='HERV')
  
  sellab <- unique(c(
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% filter(log2FoldChange>0) %>% slice_min(padj, n=5))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% filter(log2FoldChange<0) %>% slice_min(padj, n=5))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% slice_min(log2FoldChange, n=5))$loc,
    (as.data.frame(rs) %>% rownames_to_column('loc') %>% slice_max(log2FoldChange, n=5))$loc
  ))
  pointsize <- ifelse(rownames(r) %in% sellab, 2.0, 1.1)
  keyvals.shape <- ifelse(r$class=='HERV', 16, 1)
  names(keyvals.shape)[keyvals.shape==16] <- 'HERV'
  names(keyvals.shape)[keyvals.shape==1] <- 'other'
  
  EnhancedVolcano(r, 
                  lab=r$display, selectLab = sellab,
                  pointSize = pointsize,
                  drawConnectors = T,
                  col = c("grey30", "grey30", "grey30", "red2"),
                  x = 'log2FoldChange', y = 'padj',
                  labSize = 3.0,
                  shapeCustom = keyvals.shape,
                  title = n, subtitle=NULL, legendLabels=F,
                  pCutoff=p.cutoff, FCcutoff=lfc.cutoff)
})





if(makePDF) pdf(out_pdf2, paper='letter')
for(p in pl2) {
  print(p)
}
if(makePDF) dev.off()
