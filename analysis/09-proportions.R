#! /usr/bin/env Rscript

library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggpubr)

#--- Load R data
load(snakemake@input[["samp_rdata"]])
load(snakemake@input[["rtx_rdata"]])
load(snakemake@input[["metrics_rdata"]])
load(snakemake@input[["deseq_rdata"]])

if(FALSE) {
    load('analysis/01-load_sample_data.Rdata')
    load('analysis/03-load_rtx_data.Rdata')
    load('analysis/04-load_metrics.Rdata')
    load('analysis/05-deseq.Rdata')
}

#--- Sample table
samples$sample <- as.character(samples$sample_id)

#--- Transcriptome
counts.tx <- counts(dds.tx, normalized=F)
stopifnot(all(rownames(metrics) == colnames(counts.tx)))
prop.tx <- colSums(counts.tx) / metrics$ts.naln

#--- Retrotranscriptome
counts.rtx <- counts(dds.rtx, normalized=F)
stopifnot(all(rownames(metrics) == colnames(counts.rtx)))
prop.rtx <- colSums(counts.rtx) / metrics$ts.naln

#--- HERV and L1
stopifnot(all(retro.annot$locus == rownames(counts.rtx)))
prop.herv <- colSums(counts.rtx[retro.annot$class=='HERV',]) / metrics$ts.naln
prop.l1 <- colSums(counts.rtx[retro.annot$class=='L1',]) / metrics$ts.naln

#--- Data
d <- data.frame(prop.rtx=prop.rtx, 
                prop.herv=prop.herv, 
                prop.l1=prop.l1, 
                prop.tx=prop.tx
                ) %>%
    tibble::rownames_to_column('sample') %>%
    dplyr::left_join(samples, by='sample') %>%
    dplyr::mutate(prop.other = 1 - prop.tx - prop.rtx)


comps <- list( c("UVM.TP", "SKCM.TP"), c("UVM.TP", "SKCM.TM"), c("SKCM.TP", "SKCM.TM") )

pdf(snakemake@output[["pdf"]], paper='letter')

#--- RTX plot
ggboxplot(d, x="category", y="prop.rtx", fill="category") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to retrotranscriptome') +
    stat_compare_means(method="anova", label.y = log10(0.03)) +
    stat_compare_means(comparisons=comps, hide.ns=TRUE)

#--- HERV plot
ggboxplot(d, x="category", y="prop.herv", fill="category") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to HERV') +
    stat_compare_means(method="anova", label.y = log10(0.1)) +
    stat_compare_means(comparisons=comps, hide.ns=TRUE)

#--- L1 plot
ggboxplot(d, x="category", y="prop.l1", fill="category") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to L1') +
    stat_compare_means(method="anova", label.y = log10(0.01)) +
    stat_compare_means(comparisons=comps, hide.ns=TRUE)

#--- TX plot
ggboxplot(d, x="category", y="prop.tx", fill="category") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to transcriptome') +
    stat_compare_means(method="anova", label.y = log10(1.005)) +
    stat_compare_means(comparisons=comps, hide.ns=TRUE)

#--- Other plot
ggboxplot(d, x="category", y="prop.other", fill="category") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to other') +
    stat_compare_means(method="anova", label.y = log10(0.15)) +
    stat_compare_means(comparisons=comps, hide.ns=TRUE)

dev.off()
