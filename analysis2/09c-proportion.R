#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  #TODO
  # load(snakemake@input[["deseq_rdata"]])
  # load(snakemake@input[["clin_rdata"]])
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load('analysis2/03-gene_data.Rdata')  
  load("analysis2/05-filter_counts.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  
  out_pdf <- 'analysis2/09-de/proportion.pdf'
  # out_txt <- 'analysis2/09-de.txt'
}

makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)

source('analysis2/constants.R')

################################################################################
### HERV expression data
################################################################################
countDat <- mfilt.comb
cat(sprintf('%d variables\n', nrow(countDat)))

total.counts <- colSums(countDat)

byclass <- countDat %>%
  rownames_to_column("locus") %>%
  mutate(class=loc_class[locus,]$class) %>% 
  group_by(class) %>% 
  summarize(across(starts_with('TCGA'), sum))

byclass <- as.data.frame(byclass)
rownames(byclass) <- byclass$class
byclass$class <- NULL

all(colSums(byclass) == total.counts)

cpm.byclass <- t(t(byclass) / (total.counts / 1e6))
prop.byclass <- t(t(byclass) / (total.counts))

lng <- data.frame(prop.byclass) %>%
  rownames_to_column("class") %>%
  pivot_longer(!class) %>%
  mutate(name=gsub('\\.','-',name)) %>%
  mutate(clust.herv = sapply(name, function(n){clust.df[n,]$clust.retro.k4}))

pdf(out_pdf, paper='letter')

lng %>%
  filter(class == 'HERV') %>%
  ggplot(aes(x=clust.herv, y=value, color=clust.herv)) + 
  geom_boxplot(outlier.shape=NA, notch=F) + geom_jitter(position=position_jitter(0.2)) +
  scale_color_manual(values=my_pals$clust.retro.k4) +
  theme_minimal() + theme(legend.position="top") +
  labs(title="HERV proportion", x="HERV cluster", y="Proportion of reads")

lng %>%
  filter(class == 'L1') %>%
  ggplot(aes(x=clust.herv, y=value, color=clust.herv)) + 
  geom_boxplot(outlier.shape=NA, notch=F) + geom_jitter(position=position_jitter(0.2)) +
  scale_color_manual(values=my_pals$clust.retro.k4) +
  theme_minimal() + theme(legend.position="top") +
  labs(title="L1 proportion", x="HERV cluster", y="Proportion of reads")

lng %>%
  filter(class == 'protein_coding') %>%
  ggplot(aes(x=clust.herv, y=value, color=clust.herv)) + 
  geom_boxplot(outlier.shape=NA, notch=F) + geom_jitter(position=position_jitter(0.2)) +
  scale_color_manual(values=my_pals$clust.retro.k4) +
  theme_minimal() + theme(legend.position="top") +
  labs(title="Protein coding proportion", x="HERV cluster", y="Proportion of reads")

lng %>%
  filter(class == 'lncRNA') %>%
  ggplot(aes(x=clust.herv, y=value, color=clust.herv)) + 
  geom_boxplot(outlier.shape=NA, notch=F) + geom_jitter(position=position_jitter(0.2)) +
  scale_color_manual(values=my_pals$clust.retro.k4) +
  theme_minimal() + theme(legend.position="top") +
  labs(title="lncRNA proportion", x="HERV cluster", y="Proportion of reads")

dev.off()

