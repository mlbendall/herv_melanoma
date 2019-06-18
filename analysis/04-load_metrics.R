#! /usr/bin/env Rscript

library(tidyverse)

#--- Function for calculating size factors
calculateSizeFactors <- function(mapped_frags) {
    geomean <- expm1(mean(log1p(mapped_frags)))
    mapped_frags / geomean
}

#--- Load samples
load(snakemake@input[["samp_rdata"]])

cat("--------- Sample table ---------\n")
samples %>% head


#--- Load bowtie2 alignment metrics
files <- snakemake@input[["bt2_logs"]]
names(files) <- sapply(files, function(x) strsplit(x, '/')[[1]][2])
files <- files[rownames(samples)]

met.aln <- lapply(1:length(files),
               function(i){
                   lines <- readLines(files[i])
                   c(names(files)[i],
                     as.numeric(gsub('^(\\d+) .*', '\\1',  lines[1])),
                     as.numeric(gsub('^(\\d+\\.\\d+)% .*', '\\1',  lines[length(lines)]))
                   )
               }) %>% do.call(rbind, .) %>% data.frame(stringsAsFactors=F)

names(met.aln) <- c('sample', 'bt2.total_reads', 'bt2.alnrate')
met.aln$bt2.total_reads <- as.integer(met.aln$bt2.total_reads)
met.aln$bt2.alnrate <- as.double(met.aln$bt2.alnrate) * 1e-2
met.aln <- met.aln %>% mutate(bt2.naln=floor(bt2.total_reads*bt2.alnrate))
row.names(met.aln) <- met.aln$sample
met.aln$sample <- NULL

cat("-------- bowtie2 report --------\n")
met.aln %>% head

#--- Load telescope alignment metrics
files <- snakemake@input[["tele_files"]]
names(files) <- sapply(files, function(x) strsplit(x, '/')[[1]][2])
files <- files[rownames(samples)]

metrics.list <- lapply(1:length(files),
                       function(i){
                           if(file.exists(files[i])) {
                               h <- readLines(files[i], 1) %>% strsplit(., '\t') %>% unlist
                               rstr <- sapply(strsplit(h[-c(1,2)], ':'), function(t) as.numeric(unlist(t[2][1])))
                               names(rstr) <- sapply(strsplit(h[-c(1,2)], ':'), function(t) t[1])
                           } else{
                               rstr <- c(NA)
                           }
                           rstr
                       }
)
mn <- unique(do.call(c, lapply(metrics.list, names)))
met.ts <- lapply(metrics.list, function(m) {
    ret <- sapply(mn, function(x) m[x])
    names(ret) <- gsub('\\..*', '', names(ret))
    names(ret) <- paste0('ts.', names(ret))
    ret
}) %>% do.call(rbind, .) %>% data.frame
row.names(met.ts) <- rownames(samples)

cat("------- telescope report -------\n")
met.ts %>% head


#--- Combine
stopifnot(all(rownames(met.aln) == rownames(met.ts)))
stopifnot(all(met.aln$bt2.total_reads == met.ts$ts.total_fragments))

metrics <- cbind(met.aln, met.ts) %>% 
    mutate(
        total_fragments = ts.total_fragments,
        ts.naln = ts.total_fragments - ts.unmapped,
        ts.alnrate = (ts.total_fragments - ts.unmapped) / ts.total_fragments
    ) %>%
    mutate(
        bt2.size_factor = calculateSizeFactors(bt2.naln),    
        ts.size_factor = calculateSizeFactors(ts.naln)
    ) %>%
    select(total_fragments, 
           bt2.naln, bt2.alnrate, bt2.size_factor,
           ts.naln, ts.alnrate, ts.size_factor,
           ts.unique, ts.ambig, ts.overlap_unique, ts.overlap_ambig
    )

rownames(metrics) <- rownames(samples)

cat("------------ metrics -----------\n")
metrics %>% head

save(metrics, file=snakemake@output[[1]])
