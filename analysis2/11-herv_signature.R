#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(pheatmap)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
} else {
  cat("Not using snakemake\n")
  load('analysis2/03-gene_data.Rdata')
  load('analysis2/06a-dendro.Rdata')  
  load('analysis2/09-de.Rdata')  
  load("analysis2/10-feature_selection.Rdata")
  out_pdf <- 'analysis2/10-feature_selection/herv_signature.pdf'
}

source('analysis2/constants.R')

makePDF <- FALSE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)

###########################################################################################
### From 09b-heatmap.R
###########################################################################################
retro.annot$selChrom <- 'other'
retro.annot[grepl('^chr3', retro.annot$chrom), ]$selChrom <- 'chr3'
retro.annot[grepl('^chr8', retro.annot$chrom), ]$selChrom <- 'chr8'
retro.annot[grepl('^chr19', retro.annot$chrom), ]$selChrom <- 'chr19'
retro.annot$selChrom <- factor(retro.annot$selChrom, levels=c('other','chr3','chr8','chr19'))

chrom_colors <- c('other'="#efefef", 'chr3'="#0091c6", 'chr8'="#a52200", 'chr19'='#000000')


retro.annot$Group <- sapply(retro.annot$family, function(fam) retro.fam[retro.fam$family==fam,]$group)
retro.annot$Group <- factor(retro.annot$Group, levels=c('L1', names(my_pals$herv.group)))

retro.annot$selGroup <- 'other'
retro.annot[retro.annot$Group=='HERVK', ]$selGroup <- 'HERVK'
retro.annot[retro.annot$Group=='HERVH', ]$selGroup <- 'HERVH'
retro.annot[retro.annot$Group=='HERVW', ]$selGroup <- 'HERVW/9'
retro.annot[retro.annot$Group=='MER4', ]$selGroup <- 'MER4'
retro.annot[retro.annot$Group=='ERVL', ]$selGroup <- 'H/ERVL'
retro.annot[retro.annot$Group=='HERVL', ]$selGroup <- 'H/ERVL'
retro.annot[retro.annot$Group=='HERVE', ]$selGroup <- 'HERVE'
retro.annot[retro.annot$Group=='ERV3', ]$selGroup <- 'ERV3'
retro.annot[retro.annot$Group=='HUERS', ]$selGroup <- 'HUERS'

group_colors <- c(
  'other'="#efefef",
  'HERVK'="#E64B35FF",
  'HERVH'="#4DBBD5FF",
  'HERVW/9' = "#00A087FF",
  'MER4'="#3C5488FF",
  'H/ERVL'="#F39B7FFF",
  'HERVE'= "#8491B4FF" ,
  'ERV3' = "#91D1C2FF",
  'HUERS' = "#DC0000FF"
)
#### "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
retro.annot$selGroup <- factor(retro.annot$selGroup, levels=names(group_colors))

makeheatmap <- function(topgenes, ...) {
  args <- list(...)
  mat <- assay(tform)[unique(topgenes), ]
  rowdist <- as.dist((1 - cor(t(mat), method='spearman'))/2)
  
  annotation_col <- data.frame(
    row.names = colnames(dds),
    HERV.cluster = colData(dds)$clust.herv,
    stringsAsFactors=T
  )
  
  annotation_row <- data.frame(
    row.names = rownames(mat),
    Group=retro.annot[rownames(mat),]$selGroup,
    Chrom=retro.annot[rownames(mat),]$selChrom
  )
  
  annotation_colors=list(
    HERV.cluster = my_pals$clust.retro.k4,
    Group = group_colors,
    Chrom = chrom_colors
  )
  
  pheatmap(mat,
           color=viridis::inferno(20),
           scale="row", breaks=seq(-3,3,length.out=20),
           cluster_cols = as.hclust(dobj),
           clustering_distance_rows = rowdist,
           clustering_method="average",
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           annotation_row = annotation_row,
           show_colnames=F, show_rownames = T,
           fontsize = 8, fontsize_row = 6, border_color = NA,
           legend=T,
           treeheight_col=10, treeheight_row=10, cutree_cols=4,
           ...
  )
}
###########################################################################################
if(makePDF) pdf(out_pdf, paper='USr')
p <- makeheatmap(selected_vars$lasso, main="HERV signature", cellheight=9)
if(makePDF) dev.off()

library(cowplot)
library(gridExtra)

labels <- p$tree_row$labels[p$tree_row$order]

plist <- lapply(selected_vars$lasso, function(loc){
  plotCountsErrorbar(dds, loc, xaxis='clust.herv') + theme_cowplot()
})
grid.arrange(grobs=plist, nrow=3)



###########################################################################################
### What loci are signficantly different between HC3 and HC4?
###########################################################################################
subset(res$C3vC4[selected_vars$lasso,], padj<1e-3 & log2FoldChange>0)

plotCountsErrorbar(dds, 'HERVH_6q24.1b', xaxis='clust.herv')
plotCountsErrorbar(dds, 'HERVH_10q23.32', xaxis='clust.herv')
plotCountsErrorbar(dds, 'ERVLE_17p11.2c', xaxis='clust.herv')
plotCountsErrorbar(dds, 'LTR46_Xq11.1',  xaxis='clust.herv')

subset(res$C3vC4[selected_vars$lasso,], padj<1e-3 & log2FoldChange<0)

plotCountsErrorbar(dds, 'ERVLE_2p13.2a', xaxis='clust.herv')
plotCountsErrorbar(dds, 'ERVLE_5q13.2a', xaxis='clust.herv')
plotCountsErrorbar(dds, 'HERVE_Xp11.23', xaxis='clust.herv')



# split1
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] >= 10.3238 ])
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] < 10.3238])

# split2
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] >= 10.3238 ])
# table(samples$clust.herv[mat.sel['LTR46_Xq11.1',] < 10.3238])
