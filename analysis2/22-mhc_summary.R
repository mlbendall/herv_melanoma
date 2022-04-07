library(tidyverse)
library(cowplot)

orfs <- c("MER4_ORF1",
          "MER4_ORF2",
          "HML3_ORF1",
          "HML3_ORF2",
          "HML3_ORF3",
          "HML3_ORF4", 
          "HML3_ORF5",
          "HML3_ORF6",
          "HML3_ORF7",
          "HML3_ORF8",
          "HERVE_ORF1",
          "HERVE_ORF2",
          "HERVE_ORF3",
          "HERVE_ORF4",
          "HERVE_ORF5",
          "HERVE_ORF6",
          "HERVE_ORF7",
          "HERVE_ORF8",
          "HERVE_ORF9",
          "HERVE_ORF10",
          "HERVE_ORF11"
          )
gene <- sapply(strsplit(orfs, split='_'), function(x) x[1])
alleles <- c("HLA-A*02:01", # 40 
             "HLA-A*03:01", # 16
             "HLA-A*01:01", # 16
             "HLA-A*11:01", # 14
             "HLA-A*24:02", # 11
             "HLA-B*44:02", # 18
             "HLA-B*07:02", # 17
             "HLA-B*18:01", # 13
             "HLA-B*08:01"  # 13
             )


mhc_bind <- read.table('coding_analysis/mhc_binding.tsv', stringsAsFactors = F, header=T) %>%
  tidyr::separate(identity, c('locus', 'orf'), sep='_', remove=F) %>%
  mutate(bindlevel=factor(bindlevel, levels=c('.','WB','SB'))) %>%
  mutate(
    orf = factor(identity, levels=orfs),
    mhc = factor(mhc, levels=alleles)
  )

lapply(alleles, function(al) {
  sapply(orfs, function(oid){
    mhc_bind %>%
      dplyr::filter(orf == oid & mhc == al) %>%
      dplyr::filter(bindlevel == 'SB') %>% nrow
  })
}) %>% bind_cols %>% data.frame(row.names=orfs) -> strong_binders
names(strong_binders) <- alleles

View(strong_binders)

rowSums(strong_binders)

# MER4
sum(rowSums(strong_binders)[1:2])
# HML3
sum(rowSums(strong_binders)[3:10])
# HERVE
sum(rowSums(strong_binders)[11:21])


mhc_bind %>%
  dplyr::filter(bindlevel == 'SB') %>%
  select(orf, pos, mhc, peptide, score_EL, rank_EL, score_BA, rank_BA, aff) %>%
  arrange(-aff) %>% 
  write.table(file='analysis2/22-table_S2.txt', sep='\t', quote=F, row.names=F)
