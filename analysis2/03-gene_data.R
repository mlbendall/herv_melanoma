#! /usr/bin/env Rscript

library(tidyverse)

# Snakemake or not
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  ttg_tsv <- snakemake@input[["ttg_tsv"]]
  gsym_tsv <- snakemake@input[["gsym_tsv"]]
  herv_tsv <- snakemake@input[["herv_tsv"]]
  l1_tsv <- snakemake@input[["l1_tsv"]]
  tedb_tsv <- snakemake@input[['tedb_tsv']]
  geneinfo_tsv <- snakemake@input[['geneinfo_tsv']]
  out_dat <- snakemake@output[[1]]
} else {
  cat("Not using snakemake\n")  
  ttg_tsv <- 'refs/annotation/ttg.tsv'
  gsym_tsv <- 'refs/annotation/gsym.tsv'
  herv_tsv <- 'refs/annotation/HERV_rmsk.hg38.v2.tsv'
  l1_tsv <- 'refs/annotation/L1Base.hg38.v1.tsv'
  tedb_tsv <- 'refs/annotation/TE_annotation.v2.0.tsv'
  geneinfo_tsv <- 'refs/annotation/gencode.gene.info.v22.tsv'
  out_dat <- 'analysis2/03-gene_data.Rdata'
}

################################################################################
### Ensembl genes
################################################################################
#--- Load transcript to gene mapping
ttg <- read.table(ttg_tsv, 
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F
)
rownames(ttg) <- ttg$TXNAME

#--- Load gene ID to symbol mapping
gsym <- read.table(gsym_tsv, 
                   sep = '\t',
                   header = T,
                   stringsAsFactors = F
)
rownames(gsym) <- gsym$GENEID

# Use gene.info from GDC
genes.annot <- read.table(geneinfo_tsv, sep='\t', header=T)
rownames(genes.annot) <- genes.annot$gene_id
genes.annot <- genes.annot[gsym$GENEID,]


geneclasses <- recode(genes.annot$gene_type,
                      `protein_coding`="protein_coding", 
                      `miRNA`="miRNA",
                      `3prime_overlapping_ncRNA`="lncRNA",
                      `antisense`="lncRNA",
                      `bidirectional_promoter_lncRNA`="lncRNA",
                      `lincRNA`="lncRNA",
                      `macro_lncRNA`="lncRNA",
                      `non_coding`="lncRNA",
                      `processed_transcript`="lncRNA",
                      `sense_intronic`="lncRNA",
                      `sense_overlapping`="lncRNA",
                      .default="other"
)
geneclasses[grepl('pseudo', genes.annot$gene_type)] <- 'pseudo'
genes.annot$class <- factor(geneclasses, 
                          levels=c('protein_coding', 'lncRNA', 'miRNA', 'pseudo', 'other'))
table(genes.annot$class, useNA='ifany')
rm(geneclasses)

cat("---------- Tx to Gene ----------\n")
ttg %>% head
cat("--------- Gene to Sym ----------\n")
gsym %>% head
cat("--------- genes.annot ----------\n")
genes.annot %>% head



################################################################################
### Retroelement families
################################################################################
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

################################################################################
### Retroelement genes
################################################################################
annot.herv <- read.table(herv_tsv, sep='\t', header=T, stringsAsFactors=F)
annot.l1 <- read.table(l1_tsv, sep='\t', header=T, stringsAsFactors=F)
retro.annot <- rbind(annot.herv, annot.l1) %>%
  mutate(
    chrom=factor(chrom, levels=unique(annot.herv$chrom)),
    class=factor(class, levels=c('HERV', 'L1')),
    family=factor(family, levels=retro.fam$family)
  ) %>%
  dplyr::arrange(chrom, start)

row.names(retro.annot) <- retro.annot$locus

#--- Load database from Luis
te.db <- read.table(tedb_tsv,
                    sep='\t', header=T, stringsAsFactors = F)
rownames(te.db) <- te.db$Locus

cat("------- Retroelement table -----\n")
retro.annot %>% head
cat("------- TE annotation database -----\n")
te.db %>% head


################################################################################
### Full table with class
################################################################################
tmp1 <- genes.annot %>% dplyr::select(locus=gene_id, class, display=gene_name) %>% mutate(class=as.character(class))
tmp2 <- retro.annot %>% mutate(display=locus) %>% dplyr::select(locus, class, display) %>% mutate(class=as.character(class))
loc_class <- rbind(tmp1, tmp2) %>% 
  mutate(class=factor(class, levels=c('protein_coding', 'lncRNA', 'miRNA', 'pseudo', 'other', 'HERV', 'L1')))
rownames(loc_class) <- loc_class$locus

cat("----------- Full table ---------\n")
loc_class %>% head
loc_class %>% tail

save(ttg, gsym, retro.annot, retro.fam, te.db, genes.annot, loc_class, file=out_dat)
rm(list=setdiff(ls(), c('ttg', 'gsym', 'retro.annot', 'retro.fam', 'te.db', 'genes.annot', 'loc_class')))

