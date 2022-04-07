#! /usr/bin/env Rscript
require(ensembldb)

if(exists("snakemake")) {
  cat("Running with snakemake\n")  
  gtf_file <- snakemake@input[[1]]
  out_sqlite <- snakemake@output[[1]]
} else {
  cat("Not using snakemake\n")
  gtf_file <- "refs/annotation/Homo_sapiens.GRCh38.99.gtf.gz"
  out_sqlite <- "refs/annotation/Homo_sapiens.GRCh38.99.sqlite"
}

DB <- ensembldb::ensDbFromGtf(gtf_file, out_sqlite)
