#! /usr/bin/env Rscript

library(tidyverse)

tumor_lvl <- c('TP','TM')
cancer_lvl <- c('UVM','SKCM')
category_lvl <- c('UVM.TP', 'SKCM.TP', 'SKCM.TM')

samples <- read.table(snakemake@input[["samples_tsv"]],
                      sep='\t', header=T, stringsAsFactors=F) %>%
    dplyr::filter(pilot == 'True') %>% 
    dplyr::mutate(
        tumor = factor(
            dplyr::recode(sample_type,
                          `Metastatic`='TM', `Primary Tumor`='TP'
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

rownames(samples) <- samples$sample_id

cat("--------- Sample table ---------\n")
samples %>% head

save(samples, file=snakemake@output[[1]])
