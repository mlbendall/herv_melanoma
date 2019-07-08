#! /usr/bin/env Rscript

library(tidyverse)

#--- Load transcript to gene mapping
ttg <- read.table(snakemake@input[["ttg_tsv"]], 
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F
                  )

cat("---------- Tx to Gene ----------\n")
ttg %>% head

#--- Load gene ID to symbol mapping
gsym <- read.table(snakemake@input[["gsym_tsv"]], 
                   sep = '\t',
                   header = T,
                   stringsAsFactors = F
                   )

cat("--------- Gene to Sym ----------\n")
gsym %>% head


#--- Retroelement family table
retro.fam <- data.frame(
    family=c('HML1', 'HML2', 'HML3', 'HML4', 'HML5', 'HML6',
             'HERVK11', 'HERVK11D', 'HERVKC4', 'HERVK14C', 'HERVW', 'HERV9',
             'HERV30', 'HERVE', 'HERVEA', 'HERVFC1', 'HERVFC2', 'HERVFH19',
             'HERVFH21', 'HERVFRD', 'HERVH', 'HERVH48', 'HERVI', 'HERVIP10F',
             'HERVIP10FH', 'HERVL', 'HERVL18', 'HERVL32', 'HERVL40', 'HERVL66',
             'HERVL74', 'ERVL', 'ERVLB4', 'ERVLE', 'HERV3', 'HERV4',
             'HERVP71A', 'HERVS71', 'PABLA', 'PABLB', 'MER4', 'MER4B',
             'MER34B', 'MER41', 'MER61', 'MER101', 'PRIMA4', 'PRIMA41',
             'PRIMAX', 'LTR19', 'LTR23', 'LTR25', 'LTR46', 'LTR57',
             'ERV316A3', 'HARLEQUIN', 'HUERSP1', 'HUERSP2', 'HUERSP3', 'HUERSP3B',
             'L1'),
    group=c('HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVK',
            'HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVW', 'HERVW',
            'HERVW', 'HERVE', 'HERVE', 'HERVF', 'HERVF', 'HERVF',
            'HERVF', 'HERVF', 'HERVH', 'HERVH', 'HERVI', 'HERVI',
            'HERVI', 'HERVL', 'HERVL', 'HERVL', 'HERVL', 'HERVL',
            'HERVL', 'ERVL', 'ERVL', 'ERVL', 'ERV1', 'ERV1',
            'HERVP', 'HERVS', 'PAB', 'PAB', 'MER4', 'MER4',
            'MER4', 'MER4', 'MER4', 'MER4', 'PRIMA', 'PRIMA',
            'PRIMA', 'MISC', 'MISC', 'MISC', 'MISC', 'MISC',
            'ERV3', 'HARL', 'HUERS', 'HUERS', 'HUERS', 'HUERS', 
            'L1'),
    letter=c('K', 'K', 'K', 'K', 'K', 'K',
             'K', 'K', 'K', 'K', 'W', 'W',
             'W', 'E', 'E', 'F', 'F', 'F',
             'F', 'F', 'H', 'H', 'I', 'I',
             'I', 'L', 'L', 'L', 'L', 'L',
             'L', 'L', 'L', 'L', '1', '1',
             'P', 'S', '1', '1', '4', '4',
             '4', '4', '4', '4', '4', '4',
             '4', '1', '4', '4', 'H', 'L',
             '3', '4', '4', '4', '4', '4',
             'L1'),
    stringsAsFactors=F)

cat("----- Retroelement families ----\n")
retro.fam %>% head

#--- Load annotation tables
annot.herv <- read.table(snakemake@input[["herv_tsv"]], sep='\t', 
                         header=T, stringsAsFactors=F)
annot.l1 <- read.table(snakemake@input[["l1_tsv"]], sep='\t',
                       header=T, stringsAsFactors=F)

retro.annot <- rbind(annot.herv, annot.l1) %>%
    mutate(
        chrom=factor(chrom, levels=unique(annot.herv$chrom)),
        class=factor(class, levels=c('HERV', 'L1')),
        family=factor(family, levels=retro.fam$family)
    ) %>%
    dplyr::arrange(chrom, start)

row.names(retro.annot) <- retro.annot$locus

cat("------- Retroelement table -----\n")
retro.annot %>% head

save(ttg, gsym, retro.annot, retro.fam, file=snakemake@output[[1]])
