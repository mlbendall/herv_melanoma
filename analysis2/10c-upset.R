#! /usr/bin/env Rscript

library(tidyverse)
library(UpSetR)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
} else {
  cat("Not using snakemake\n")
  load("analysis2/10-feature_selection.Rdata")
  out_pdf <- 'analysis2/10-feature_selection/upset.pdf'
}


dir.create(dirname(out_pdf), showWarnings = F)
pdf(out_pdf, width=6, height=3)
upset(fromList(selected_vars), order.by = "freq")
dev.off()

