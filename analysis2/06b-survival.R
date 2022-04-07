#! /usr/bin/env Rscript

library(tidyverse)
library(survival)
library(survminer)


#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["clin_rdata"]])
  load(snakemake@input[["clust_rdata"]])
  out_pdf <- snakemake@output[["outpdf"]]
  
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  out_pdf <- 'analysis2/06-unsupervised/survival.pdf'
}

# library(TCGAbiolinks)
# clin.uvm <- GDCquery_clinic("TCGA-UVM", "clinical")
# rownames(clin.uvm) <- clin.uvm$submitter_id
# clin.uvm <- clin.uvm[samples$case_id,]

makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)
if(makePDF) pdf(out_pdf, paper='letter', onefile=T)
source('analysis2/constants.R')

################################################################################
### Fit 1: HERV
################################################################################
uvmdat <- mdata[,c("death", "cause_of_death", "time_to_death", "time_to_met", "time_to_met_or_death", 
                   "metastatic", "chr3CN", "D3M3", "clust.mRNA", "clust.methyl", "clust.lncRNA", "clust.miRNA")]

stopifnot(all(rownames(uvmdat)==rownames(clust.df)))
uvmdat$clust.herv <- clust.df$clust.retro.k4
uvmdat$group.herv <- clust.df$group.k4

uvmdat$time <- uvmdat$time_to_death
uvmdat$status <- ifelse(uvmdat$death %in% c("A", "AM"), 0, 1)

# Filter uvmdat
uvmdat <- uvmdat[uvmdat$death %in% c('A','AM','DM', 'DU'),]

surv_object <- Surv(uvmdat$time, uvmdat$status)
fit1 <- survfit(surv_object ~ clust.herv, data = uvmdat)

surv_pvalue(fit1)

ggsurvplot(fit1, data=uvmdat, pval=T,
           risk.table=F, ncensor.plot=F,
           xscale='d_m', break.x.by = 365.25,
           xlim=c(0, 365.25*5),
           xlab="Time (months)",
           ggtheme = theme_light(),
           conf.int = F, conf.int.style = "step",
           legend.labs=names(my_pals$clust.retro.k4),
           palette=unname(my_pals$clust.retro.k4)
)

# fit1.cox <- coxph(surv_object ~ group.herv, data = uvmdat, ties='exact')
# ggforest(fit1.cox, data=uvmdat)

rm(uvmdat, surv_object)


################################################################################
### Fit 2: D3M3
################################################################################
uvmdat <- mdata[,c("death", "cause_of_death", "time_to_death", "time_to_met", "time_to_met_or_death", 
                   "metastatic", "chr3CN", "D3M3", "clust.mRNA", "clust.methyl", "clust.lncRNA", "clust.miRNA")]

stopifnot(all(rownames(uvmdat)==rownames(clust.df)))
uvmdat$clust.herv <- clust.df$clust.retro.k4

uvmdat$time <- uvmdat$time_to_death
uvmdat$status <- ifelse(uvmdat$death %in% c("A", "AM"), 0, 1)

# Filter uvmdat
uvmdat <- uvmdat[uvmdat$death %in% c('A','AM','DM', 'DU') & uvmdat$chr3CN %in% c('1','2'), ]
uvmdat <- droplevels(uvmdat)

surv_object <- Surv(uvmdat$time, uvmdat$status)
fit2 <- survfit(surv_object ~ D3M3, data = uvmdat)
ggsurvplot(fit2, data=uvmdat, pval=T,
           risk.table=F, ncensor.plot=F,
           xscale='d_m', break.x.by = 365.25,
           xlim=c(0, 365.25*5),
           xlab="Time (months)",
           ggtheme = theme_light(),
           conf.int = F, conf.int.style = "step",
           legend.labs=names(my_pals$D3M3)[2:3],
           palette=unname(my_pals$D3M3)[2:3]
)

rm(uvmdat, surv_object)

if(makePDF) dev.off()

