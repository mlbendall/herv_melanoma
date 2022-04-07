#! /usr/bin/env Rscript

library(tidyverse)
library(ensembldb)

# Snakemake or not
if(exists("snakemake")) {
    cat("Running with snakemake\n")
    ttg_tsv <- snakemake@input[["ttg_tsv"]]
    gsym_tsv <- snakemake@input[["gsym_tsv"]]
    herv_tsv <- snakemake@input[["herv_tsv"]]
    l1_tsv <- snakemake@input[["l1_tsv"]]
    tedb_tsv <- snakemake@input[['tedb_tsv']]
    ens_sqlite <- snakemake@input[['ens_sqlite']]
    out_rdata <- snakemake@output[[1]]
} else {
    ttg_tsv <- 'refs/annotation/ttg.tsv'
    gsym_tsv <- 'refs/annotation/gsym.tsv'
    herv_tsv <- 'refs/annotation/HERV_rmsk.hg38.v2.tsv'
    l1_tsv <- 'refs/annotation/L1Base.hg38.v1.tsv'
    tedb_tsv <- 'refs/annotation/TE_annotation.v2.0.tsv'
    ens_sqlite <- 'refs/annotation/Homo_sapiens.GRCh38.99.sqlite'
    out_rdata <- 'analysisUM/01-load_gene_data.Rdata'
}

#--- Load transcript to gene mapping
ttg <- read.table(ttg_tsv, 
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F
                  )
rownames(ttg) <- ttg$TXNAME


cat("---------- Tx to Gene ----------\n")
ttg %>% head

#--- Load gene ID to symbol mapping
gsym <- read.table(gsym_tsv, 
                   sep = '\t',
                   header = T,
                   stringsAsFactors = F
                   )
rownames(gsym) <- gsym$GENEID

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
annot.herv <- read.table(herv_tsv, sep='\t', 
                         header=T, stringsAsFactors=F)
annot.l1 <- read.table(l1_tsv, sep='\t',
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


cat("------- TE annotation database -----\n")
te.db <- read.table(tedb_tsv,
                   sep='\t', header=T, stringsAsFactors = F)
rownames(te.db) <- te.db$Locus
te.db %>% head



# build_ens_DB<-function(){
#     system("curl -q ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz --output refs/annotation/Homo_sapiens.GRCh38.99.gtf.gz")
#     DB<-ensembldb::ensDbFromGtf("refs/annotation/Homo_sapiens.GRCh38.99.gtf.gz")
# }
# ifelse(!file.exists("refs/annotation/Homo_sapiens.GRCh38.99.sqlite"),
#        yes=build_ens_DB(), 
#        no= "Annotations already downloaded")

hg38_ens99<- ensembldb::EnsDb(ens_sqlite)

ucsc_style <- mapSeqlevels(seqlevels(genes(hg38_ens99,filter=SeqNameFilter(c(1:22,"X","Y")))), "UCSC")
hg38_genes <- renameSeqlevels(genes(hg38_ens99,filter=SeqNameFilter(c(1:22,"X","Y")), columns=c("gene_name","gene_biotype")), ucsc_style)
# hg38_exons <- renameSeqlevels(exons(hg38_ens99, columns = c("tx_id", "tx_biotype", "gene_name", "gene_id"),filter =SeqNameFilter(c(1:22,"X","Y"))), ucsc_style)

# hg38_5utr <- renameSeqlevels(fiveUTRsByTranscript(hg38_ens99, columns = c("tx_id", "tx_biotype", "gene_name", "gene_id"),filter =SeqNameFilter(c(1:22,"X","Y"))), ucsc_style)
# hg38_3utr <- renameSeqlevels(threeUTRsByTranscript(hg38_ens99, columns = c("tx_id", "tx_biotype", "gene_name", "gene_id"),filter =SeqNameFilter(c(1:22,"X","Y"))), ucsc_style)

save(ttg, gsym, retro.annot, retro.fam, te.db, hg38_genes, file=out_rdata)
