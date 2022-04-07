#! /usr/bin/env Rscript

library(tidyverse)

if(exists("snakemake")) {
  cat("Running with snakemake\n")  
  sample_tsv <- snakemake@input[["samples_tsv"]]
  out_dat <- snakemake@output[[1]]
} else {
  cat("Not using snakemake\n")  
  sample_tsv <- 'metadata/tcga_samples.tsv'
  out_dat <- 'analysis2/01-sample_data.Rdata'
}

tumor_lvl <- c('TP','TM', 'NT', 'TAM')
cancer_lvl <- c('UVM','SKCM')
category_lvl <- c('UVM.TP', 'SKCM.TP', 'SKCM.TM', 'SKCM.NT', 'SKCM.TAM')

samples <- read.table(sample_tsv,
                      sep='\t', header=T, stringsAsFactors=F) %>%
  dplyr::mutate(
    tumor = factor(
      dplyr::recode(sample_type,
                    `Metastatic`='TM', `Primary Tumor`='TP',
                    `Solid Tissue Normal`='NT', `Additional Metastatic`='TAM'
      ),
      levels=tumor_lvl
    )
  ) %>%
  tidyr::separate(project_id, c('x', 'cancer'), sep='-', remove=F) %>% 
  dplyr::mutate(cancer=factor(cancer, levels=cancer_lvl)) %>%
  dplyr::mutate(category=paste(cancer,tumor,sep='.')) %>%
  dplyr::mutate(category=factor(category, levels=category_lvl)) %>%
  dplyr::arrange(category) %>%
  dplyr::mutate(
    case_id=factor(case_id),
    sample_id=factor(sample_id)
  ) %>%
  dplyr::select(cancer, tumor, category, case_id, sample_id)

# Filter only UM
samples <- samples %>% filter(cancer=='UVM') %>% droplevels
rownames(samples) <- samples$sample_id

cat("--------- Sample table ---------\n")
samples %>% head

save(samples, file=out_dat)
rm(list=setdiff(ls(), c('samples')))