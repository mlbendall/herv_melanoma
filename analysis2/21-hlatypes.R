#! /usr/bin/env Rscript

library(tidyverse)

if(exists("snakemake")) {
  cat("Running with snakemake\n")  
} else {
  cat("Not using snakemake\n")  
  load('analysis2/01-sample_data.Rdata')
  load('analysis2/06-unsupervised.Rdata')
  hla_csv <- 'metadata/OptiTypeCallsHLA_20171207.tsv'
}

source('analysis2/constants.R')

hla <- read.csv(hla_csv) %>%
  tidyr::separate(aliquot_id, remove=F, c('i1','i2','i3','i4'), extra='drop') %>%
  tidyr::unite(case_id, i1:i3, sep='-', remove=F) %>% 
  tidyr::unite(sample_id, i1:i4, sep='-')


hla.uvm <- left_join(data.frame(case_id=samples$case_id), hla)
rownames(hla.uvm) <- rownames(samples)
hla.uvm <- hla.uvm[,c('A1','A2','B1','B2','C1','C2')]

alleles.A <- sort(unique(c(hla.uvm$A1, hla.uvm$A2)))
one_hot.A <- lapply(alleles.A, function(x) {
  apply(hla.uvm, 1, function(r) ifelse(r[1] == x | r[2] == x, 1, 0))
}) %>% bind_cols %>% data.frame
names(one_hot.A) <- alleles.A
rownames(one_hot.A) <- rownames(hla.uvm)

alleles.B <- sort(unique(c(hla.uvm$B1, hla.uvm$B2)))
one_hot.B <- lapply(alleles.B, function(x) {
  apply(hla.uvm, 1, function(r) ifelse(r[3] == x | r[4] == x, 1, 0))
}) %>% bind_cols %>% data.frame
names(one_hot.B) <- alleles.B
rownames(one_hot.B) <- rownames(hla.uvm)

alleles.C <- sort(unique(c(hla.uvm$C1, hla.uvm$C2)))
one_hot.C <- lapply(alleles.C, function(x) {
  apply(hla.uvm, 1, function(r) ifelse(r[5] == x | r[6] == x, 1, 0))
}) %>% bind_cols %>% data.frame
names(one_hot.C) <- alleles.C
rownames(one_hot.C) <- rownames(hla.uvm)

one_hot.uvm <- cbind(one_hot.A, one_hot.B, one_hot.C)

hla.uvm$A1 <- factor(hla.uvm$A1)
hla.uvm$A2 <- factor(hla.uvm$A2)
hla.uvm$B1 <- factor(hla.uvm$B1)
hla.uvm$B2 <- factor(hla.uvm$B2)
hla.uvm$C1 <- factor(hla.uvm$C1)
hla.uvm$C2 <- factor(hla.uvm$C2)


tmppca <- pca.obj
stopifnot(all(rownames(tmppca$metadata) == rownames(hla.uvm)))
tmppca$metadata <- cbind(tmppca$metadata, hla.uvm)
tmppca$metadata <- cbind(tmppca$metadata, lapply(one_hot.uvm, as.logical))


PCAtools::biplot(tmppca,
                 colby = 'A1',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-A 1")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'A2',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-A 2")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'B1',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-B 1")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'B2',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-B 2")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'C1',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-C 1")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'C2',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-C 2")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)


PCAtools::eigencorplot(tmppca,
                       components = PCAtools::getComponents(pca.obj, 1:4),
                       metavars = alleles.A,
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
PCAtools::biplot(tmppca,
                 colby = 'A*31:01',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-A*31:01")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)

PCAtools::biplot(tmppca,
                 colby = 'A*01:01',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-A*01:01")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)



PCAtools::eigencorplot(tmppca,
                       components = PCAtools::getComponents(pca.obj, 1:4),
                       metavars = alleles.B,
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
PCAtools::biplot(tmppca,
                 colby = 'B*44:02',
                 ellipse = FALSE,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HLA-B*44:02")
) + xlim(-45, 45) + ylim(-30, 30) + theme(aspect.ratio=1)



PCAtools::eigencorplot(tmppca,
                       components = PCAtools::getComponents(pca.obj, 1:4),
                       metavars = alleles.C,
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)



colSums(one_hot.A, na.rm=T) %>% sort
colSums(one_hot.B, na.rm=T) %>% sort


