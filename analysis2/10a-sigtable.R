#! /usr/bin/env Rscript

library(tidyverse)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
  # load(snakemake@input[["deseq_rdata"]])
  # load(snakemake@input[["clin_rdata"]])
} else {
  cat("Not using snakemake\n")
  load("analysis2/10-feature_selection.Rdata")
  load("analysis2/03-gene_data.Rdata")
}



tmp <- te.db[selected_vars$lasso,] %>%
  select(Locus, TE_CODING, TE_type, IntersectedGene, ClosestUpstream_gn, DistUpstream, ClosestDownstream_gn, DistDownstream )

tmp$Gene <- paste0(tmp$IntersectedGene, '*')
tmp[tmp$IntersectedGene=='None',]$Gene <- ifelse(
  tmp[tmp$IntersectedGene=='None',]$DistUpstream < tmp[tmp$IntersectedGene=='None',]$DistDownstream, tmp[tmp$IntersectedGene=='None',]$ClosestUpstream_gn, tmp[tmp$IntersectedGene=='None',]$ClosestDownstream_gn
  )

tmp %>% 
  mutate(region=tolower(TE_type),
         coding=recode(TE_CODING, `coding`='Y', `noncoding`='N')) %>%
  select(Locus, Gene, region,  coding) %>%
  write.table(sep='\t', quote=F, row.names = F)
