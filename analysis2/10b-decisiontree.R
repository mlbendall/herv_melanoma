#! /usr/bin/env Rscript

library(tidyverse)
library(rpart)
library(rpart.plot)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
  # load(snakemake@input[["deseq_rdata"]])
  # load(snakemake@input[["clin_rdata"]])
} else {
  cat("Not using snakemake\n")
  load("analysis2/10-feature_selection.Rdata")
  out_pdf <- 'analysis2/10-feature_selection/decisiontree.pdf'
}

source('analysis2/constants.R')

makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)

mat.sel <- mat[selected_vars$lasso,]
rpart.fit <- rpart(y ~ . , data=data.frame(y=samples$clust.herv, t(mat.sel)), method = 'class')

if(makePDF) pdf(out_pdf, paper='USr')
rpart.plot(rpart.fit, type=5, extra=1, digits=3, box.palette = as.list(my_pals$clust.retro.k4))
# rpart.plot(rpart.fit, type=5, extra=4)
if(makePDF) dev.off()



# split1
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] >= 10.3238 ])
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] < 10.3238])

# split2
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] >= 10.3238 ])
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] < 10.3238])
