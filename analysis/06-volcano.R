#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)

library(EnhancedVolcano)

#--- Parameters
padj_cutoff <- 0.1  # Adjusted p-value (FDR) cutoff
lfc_cutoff <- 1.0   # log2FoldChange (LFC) cutoff

#--- Load DESeq objects
load(snakemake@input[["deseq_rdata"]])

#--- Load gene data
load(snakemake@input[["gene_rdata"]])

cat("-------- Gene to Symbol --------\n")
gsym %>% head

pdf(snakemake@output[[1]], paper='letter')

for(i in 1:3) {
    res_sym <- results.tx[[i]] %>% 
        data.frame %>%
        tibble::rownames_to_column('row') %>%
        dplyr::inner_join(gsym, by=c("row"="GENEID")) %>%
        dplyr::select(row, SYM)

    print(
        EnhancedVolcano(results.tx[[i]],
            lab = res_sym$SYM,
            x = 'log2FoldChange',
            y = 'padj',
            xlim = c(-8, 8),
            title = paste(names(results.tx)[i], '(transcriptome)'),
            pCutoff = padj_cutoff,
            FCcutoff = lfc_cutoff,
            transcriptPointSize = 1.5,
            transcriptLabSize = 3.0
        )
    )
}

for(i in 1:3) {
    print(
        EnhancedVolcano(results.rtx[[i]],
            lab = rownames(results.rtx[[i]]),
            x = 'log2FoldChange',
            y = 'padj',
            xlim = c(-8, 8),
            title = paste(names(results.rtx)[i], '(retrotranscriptome)'),
            pCutoff = padj_cutoff,
            FCcutoff = lfc_cutoff,
            transcriptPointSize = 1.5,
            transcriptLabSize = 3.0
        )
    )
}

dev.off()
