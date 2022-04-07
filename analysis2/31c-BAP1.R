library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(cowplot)

load("analysis2/01-sample_data.Rdata")  
load("analysis2/02-clinical_data.Rdata")
load('analysis2/03-gene_data.Rdata')  
load("analysis2/05-filter_counts.Rdata")
load("analysis2/06-unsupervised.Rdata")

outdir <- 'analysis2/31-additional_counts'  
makePDF <- TRUE
if(makePDF) dir.create(outdir, showWarnings = F)

source('analysis2/constants.R')

cutoff.count <- 5
p.cutoff <- 1e-3
lfc.cutoff <- 1.5

add_labels <- c(
  "HERVH_1p31.3d", "ERV316A3_6p21.33c", "MER4_22q12.3", "MER4B_8q21.11", "HERVH_5p15.33",
  "HERVL_21q21.3f", "HML2_10q24.2", "HML2_7p22.1", "HML6_19q13.41e", "MER4B_3q21.2", "ERV316A3_3q13.12c"
)

################################################################################
### TCGA BAP1 comparison
################################################################################
countDat <- mfilt.herv
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(samples)))
stopifnot(all(rownames(mdata) == rownames(samples)))
stopifnot(all(rownames(samples) == rownames(clust.df)))

samples$mut.BAP1 <- factor(ifelse(mdata$mut.BAP1==1, 'mut','wt'), levels=c('wt','mut'))
samples$clust.herv <- clust.df$clust.retro.k4

dds.tcga <- DESeq2::DESeqDataSetFromMatrix(countDat, samples, ~mut.BAP1)
dds.tcga <- DESeq2::DESeq(dds.tcga, parallel=TRUE)
res.tcga <- results(dds.tcga, contrast=c('mut.BAP1', 'mut','wt'), alpha=p.cutoff)
res.tcga$class <- loc_class[rownames(res.tcga),]$class
res.tcga$display <- loc_class[rownames(res.tcga),]$display
sig.tcga <- subset(res.tcga, padj<p.cutoff & abs(log2FoldChange) > lfc.cutoff)


sellab1 <- unique(c(
  (as.data.frame(sig.tcga) %>% rownames_to_column('loc') %>% filter(log2FoldChange>0) %>% slice_min(padj, n=5))$loc,
  (as.data.frame(sig.tcga) %>% rownames_to_column('loc') %>% filter(log2FoldChange<0) %>% slice_min(padj, n=5))$loc,
  (as.data.frame(sig.tcga) %>% rownames_to_column('loc') %>% slice_min(log2FoldChange, n=5))$loc,
  (as.data.frame(sig.tcga) %>% rownames_to_column('loc') %>% slice_max(log2FoldChange, n=5))$loc
))
sellab1 <- unique(c(sellab1, add_labels))

g1 <- EnhancedVolcano(res.tcga, 
                lab=res.tcga$display,
                selectLab = sellab1,
                xlim=c(-6.5,6.5),
                ylim=c(0,40),
                pointSize = 1.1,
                drawConnectors = F,
                widthConnectors = 0.5,
                col = c("grey30", "grey30", "grey30", "red2"),
                x = 'log2FoldChange', y = 'padj',
                labSize = 3.5,
                axisLabSize = 8,
                title = 'TCGA_UVM', subtitle=NULL, 
                legendPosition='none',
                xlab=NULL, ylab=NULL, caption = NULL,
                pCutoff=p.cutoff, FCcutoff=lfc.cutoff) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

if(makePDF) pdf(file.path(outdir, 'tcga_bap1.pdf'), paper='letter')
print(g1)
if(makePDF) dev.off()



################################################################################
### GSE149920 BAP1 comparison
################################################################################

load("analysis2/31-additional_counts.Rdata")

stopifnot(all(rownames(samples.add) == names(counts.add.comb)))

study_id <- 'GSE149920'
in_study <- samples.add$study==study_id

samp.study <- samples.add[in_study,]

counts.study.herv <- counts.add.herv[,in_study]
mfilt.study.herv <- counts.study.herv[rowSums(counts.study.herv) > cutoff.count, ]

dds.study <- DESeq2::DESeqDataSetFromMatrix(mfilt.study.herv, samp.study, ~category)
dds.study <- DESeq2::DESeq(dds.study, parallel=T)

res.study <- results(dds.study, contrast=c('category','BAP1mut', 'BAP1wt'), alpha=p.cutoff)
res.study$class <- loc_class[rownames(res.study),]$class
res.study$display <- loc_class[rownames(res.study),]$display
sig.study <- subset(res.study, padj<p.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig.study[order(sig.study$padj),]

sellab2 <- unique(c(
  (as.data.frame(sig.study) %>% rownames_to_column('loc') %>% filter(log2FoldChange>0) %>% slice_min(padj, n=5))$loc,
  (as.data.frame(sig.study) %>% rownames_to_column('loc') %>% filter(log2FoldChange<0) %>% slice_min(padj, n=5))$loc,
  (as.data.frame(sig.study) %>% rownames_to_column('loc') %>% slice_min(log2FoldChange, n=5))$loc,
  (as.data.frame(sig.study) %>% rownames_to_column('loc') %>% slice_max(log2FoldChange, n=5))$loc
))
sellab2 <- unique(c(sellab2, add_labels))

g2 <- EnhancedVolcano(res.study, 
                lab=res.study$display,
                selectLab = sellab2,
                xlim=c(-26, 26),
                ylim=c(0,40),
                pointSize = 1.1,
                drawConnectors = F,
                widthConnectors = 0.5,
                col = c("grey30", "grey30", "grey30", "red2"),
                x = 'log2FoldChange', y = 'padj',
                labSize = 3.5,
                axisLabSize = 8,
                title = study_id, subtitle=NULL, 
                legendPosition='none',
                xlab=NULL, ylab=NULL, caption = NULL,
                pCutoff=p.cutoff, FCcutoff=lfc.cutoff) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

if(makePDF) pdf(file.path(outdir, 'GSE149920_bap1.pdf'), paper='letter')
print(g2)
if(makePDF) dev.off()


################################################################################
### compare datasets
################################################################################
load('analysis2/10-feature_selection-selected_vars.Rdata')
table(samples$mut.BAP1, samples$clust.herv)

tcga.up <- as.data.frame(sig.tcga) %>%
  rownames_to_column('loc') %>% 
  filter(log2FoldChange>0) %>%
  arrange(-log2FoldChange)

tcga.dn <- as.data.frame(sig.tcga) %>%
  rownames_to_column('loc') %>% 
  filter(log2FoldChange<0) %>%
  arrange(log2FoldChange)


study.up <- as.data.frame(sig.study) %>%
  rownames_to_column('loc') %>% 
  filter(log2FoldChange>0) %>%
  arrange(-log2FoldChange)

study.dn <- as.data.frame(sig.study) %>%
  rownames_to_column('loc') %>% 
  filter(log2FoldChange<0) %>%
  arrange(log2FoldChange)

intr.up <- intersect(tcga.up$loc, study.up$loc)
intr.dn <- intersect(tcga.dn$loc, study.dn$loc)

plotCountsErrorbar(dds.tcga, "MER4_22q12.3" , 'mut.BAP1', color = 'clust.herv') + theme_cowplot()
plotCountsErrorbar(dds.study, "HML6_19q13.41e", 'case_id') + theme_cowplot()



plotCountsErrorbar(dds.tcga, "MER4B_3q21.2", 'mut.BAP1', color = 'clust.herv') + theme_cowplot()
plotCountsErrorbar(dds.tcga, "ERV316A3_6p21.33c", 'mut.BAP1', color = 'clust.herv') + theme_cowplot()


plotCountsErrorbar(dds.tcga, 'HERVH_17q25.1c', 'mut.BAP1', color = 'clust.herv') + theme_cowplot()


  [sig.tcga$log2FoldChange>0, ]


plotCountsErrorbar(dds.tcga, 'ERV316A3_3q13.12c', 'case_id')



# plotCountsErrorbar(dds.study, 'ERV316A3_3q13.12c', 'case_id')
# plotCountsErrorbar(dds.study, 'MER4B_3q21.2', 'case_id')