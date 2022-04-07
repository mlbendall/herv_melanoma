#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(PCAtools)
library(RColorBrewer)
library(pheatmap)
library(TCGAbiolinks)

#--- Load DESeq objects
if(exists("snakemake")) {
    cat("Running with snakemake\n")
    load(snakemake@input[["deseq_rdata"]])
    load(snakemake@input[["clin_rdata"]])
} else {
    cat("Not using snakemake\n")    
    load("analysisUM/05-deseq_int.Rdata")
    load("analysisUM/01-load_clin_data.Rdata")
}

stopifnot(all(colData(tform.rtx)$sample_id == colData(tform.tx)$sample_id))
stopifnot(all(colData(tform.rtx)$sample_id == row.names(mdata)))

if(exists("snakemake")) {
    pdf(snakemake@output[[1]], paper='letter')
} else {
    pdf("analysisUM/10-allsamples-pca.pdf", paper='letter')
}

#--- Retrotranscriptome
p.rtx <- PCAtools::pca(assay(tform.rtx), metadata=mdata, removeVar = 0.50)





### Screeplot
varline <- 50
varline.x <- min(which(cumsum(p.rtx$variance) >= varline))

PCAtools::screeplot(p.rtx,
    components = getComponents(p.rtx, 1:30),
    title="Retrotranscriptome SCREE",
    hline=varline, vline=varline.x
) +
geom_text(aes(5, varline, label = paste0(varline, '% explained variation'), vjust = -1))

### Loadings
rangeRetain <- 0.01
PCAtools::plotloadings(p.rtx,
    title=paste0("Retrotranscriptome loadings"),
    rangeRetain = rangeRetain,
    caption = paste0('Top ', rangeRetain * 100, '% variables'),
    subtitle = 'PC1-PC5',
    shapeSizeRange = c(3,3),    
    labSize = 3.0,    
    shape = 24,
    col = c('limegreen', 'black', 'red3'),
    drawConnectors = TRUE
)

### PCA plot wrappers
biplot_wrapper <- function(color.by, colkey=NULL) {
    PCAtools::biplot(p.rtx,
        colby = color.by,
        colkey = colkey,    
        hline = 0,
        vline = 0,
        legendPosition = 'right',
        lab=F,
        title=paste0("Retrotranscriptome PCA (", color.by, ")")
    )
}

pairplot_wrapper <- function(maxcomp, color.by, colkey) {
    PCAtools::pairsplot(p.rtx, 
        components = getComponents(p.rtx, c(1:maxcomp)),
        colby = color.by,
        colkey = colkey,
        pointSize = 0.4,
        title=paste0("Retrotranscriptome PC1-PC", maxcomp, " (", color.by, ")")
    )
}

#SCNA cluster coloring
biplot_wrapper('clust.SCNA', c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'))
pairplot_wrapper(5, 'clust.SCNA', c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'))

# Metastatic colorings
biplot_wrapper('metastatic', c("No"="#666666", "Yes"="#990000"))
pairplot_wrapper(5, 'metastatic', c("No"="#666666", "Yes"="#990000"))

# BAP1 coloring
biplot_wrapper("BAP1", c("NA"="#CCCCCC", "NON"="#a60e8b", "MIS"="#ff77ad", "FS"="#df0192", "SPL"="#8981c0", "IDEL"="#ffbec3", "HDEL"="#ffefff"))
pairplot_wrapper(5, "BAP1", c("NA"="#CCCCCC", "NON"="#a60e8b", "MIS"="#ff77ad", "FS"="#df0192", "SPL"="#8981c0", "IDEL"="#ffbec3", "HDEL"="#ffefff"))

# GNAQ coloring
biplot_wrapper("GNAQ", c("NA"="#CCCCCC", "G48L"="#a60e8b", "Q209L"="#ff77ad", "Q209P"="#df0192", "R183Q"="#8981c0"))
pairplot_wrapper(5, "GNAQ", c("NA"="#CCCCCC", "G48L"="#a60e8b", "Q209L"="#ff77ad", "Q209P"="#df0192", "R183Q"="#8981c0"))

# GNA11 coloring
biplot_wrapper("GNA11", c("NA"="#CCCCCC", "Q209L"="#ff77ad", "R183C"="#8981c0", "R166H"="#9981c0"))
pairplot_wrapper(5, "GNA11", c("NA"="#CCCCCC", "Q209L"="#ff77ad", "R183C"="#8981c0", "R166H"="#9981c0"))


# Metastatic with labels
color.by <- 'metastatic'
colkey <- c("No"="#666666", "Yes"="#990000")
PCAtools::biplot(p.rtx,
    colby = color.by,
    colkey = colkey,
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=T,
    selectLab = p.rtx$yvars[p.rtx$rotated$PC2 < 0 & mdata$metastatic=='Yes'],
    title=paste0("Retrotranscriptome PCA (", color.by, ")")
)

p.rtx$metadata$clust.rtxPC1 <- factor(ifelse(p.rtx$rotated$PC1 > 0, 'PC1a', 'PC1b'), levels=c('PC1a','PC1b'))
p.rtx$metadata$clust.rtxPC2 <- factor(ifelse(p.rtx$rotated$PC2 > 0, 'PC2a', 'PC2b'), levels=c('PC2a','PC2b'))
p.rtx$metadata$clust.rtxPC3 <- factor(ifelse(p.rtx$rotated$PC3 > -40, 'PC3a', 'PC3b'), levels=c('PC3a','PC3b'))

color.by <- "clust.rtxPC1"
PCAtools::biplot(p.rtx,
    colby = color.by,
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=T,
    selectLab = p.rtx$yvars[p.rtx$rotated$PC2 < 0 & mdata$metastatic=='Yes'],
    title=paste0("Retrotranscriptome PCA (", color.by, ")")
)

color.by <- "clust.rtxPC2"
PCAtools::biplot(p.rtx,
    colby = color.by,
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=T,
    selectLab = p.rtx$yvars[p.rtx$rotated$PC2 < 0 & mdata$metastatic=='Yes'],
    title=paste0("Retrotranscriptome PCA (", color.by, ")")
)

color.by <- "clust.rtxPC3"
PCAtools::biplot(p.rtx,
    y = 'PC3',
    colby = color.by,
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=T,
    selectLab = p.rtx$yvars[p.rtx$rotated$PC2 < 0 & mdata$metastatic=='Yes'],
    title=paste0("Retrotranscriptome PCA (", color.by, ")")
)

dev.off()

# DE genes based on PC2

# Recalculate DESeq object
stopifnot(rownames(colData(dds.rtx)) == p.rtx$yvars)
colData(dds.rtx)$clust.rtxPC1 <- factor(ifelse(p.rtx$rotated$PC1 > 0, 'PC1a', 'PC1b'), levels=c('PC1a','PC1b'))
colData(dds.rtx)$clust.rtxPC2 <- factor(ifelse(p.rtx$rotated$PC2 > 0, 'PC2a', 'PC2b'), levels=c('PC2a','PC2b'))
colData(dds.rtx)$clust.rtxPC3 <- factor(ifelse(p.rtx$rotated$PC3 > -40, 'PC3a', 'PC3b'), levels=c('PC3a','PC3b'))

design(dds.rtx) <- ~clust.rtxPC1+clust.rtxPC2
dds.rtx <- DESeq(dds.rtx, parallel=T)

# Criteria
alpha <- 0.01
lfc <- 1
ntop.hm <- 50
summary(res1 <- results(dds.rtx, contrast=c("clust.rtxPC1", "PC1a", "PC1b"), alpha=alpha))
summary(res2 <- results(dds.rtx, contrast=c("clust.rtxPC2", "PC2a", "PC2b"), alpha=alpha))

sig2 <- subset(res2, padj<alpha & abs(log2FoldChange)>=lfc)
nameorder <- rownames(sig2[order(sig2$padj),])

# gets only HERV
topgenes <- nameorder[grep('^L1', nameorder,invert=T)][1:ntop.hm]
mat <- assay(tform.rtx)[topgenes,]

stopifnot(all(rownames(mdata) == colnames(mat)))
rowdist <- as.dist((1 - cor(t(mat), method='spearman'))/2)

### HEATMAP

if(exists("snakemake")) {
    pdf(snakemake@output[[2]], paper='letter')
} else {
    pdf("analysisUM/07-pca_heatmap.pdf", paper='letter')
}


renamed_clusters <- plyr::revalue(colData(dds.rtx)$clust.rtxPC2, c("PC2a"="A", "PC2b"="B"))

annotation_col <- data.frame(
    row.names = rownames(mdata),
    clust.methyl=mdata$clust.methyl,
    clust.mRNA=mdata$clust.mRNA,
    clust.SCNA=mdata$clust.SCNA,
    BAP1=mdata$BAP1,
    metastatic=mdata$metastatic,
    clust.HERV = renamed_clusters,
    stringsAsFactors=T
)


cA <- brewer.pal(n = 3, name = "RdBu")[3]
cB <- brewer.pal(n = 3, name = "RdBu")[1]

pheatmap(mat,
         scale="row", 
         clustering_distance_rows = rowdist,
         clustering_method="average",
         annotation_col = annotation_col,
         annotation_colors=list(
             clust.HERV = c("A"=cA, "B"=cB),
             metastatic = c("No"="#666666", "Yes"="#990000"),
             BAP1 = c("NA"="#CCCCCC", "NON"="#a60e8b", "MIS"="#ff77ad", "FS"="#df0192", "SPL"="#8981c0", "IDEL"="#ffbec3", "HDEL"="#ffefff"),
             clust.SCNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
             clust.mRNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
             clust.methyl = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6')
         ),
         show_colnames=F,
         cutree_cols = 2,
         fontsize = 8, fontsize_row = 6, border_color = NA
)

dev.off()


clin.uvm <- GDCquery_clinic("TCGA-UVM", "clinical")

surv_clust <- renamed_clusters # colData(dds.rtx)$clust.rtxPC2
names(surv_clust) <- sapply(strsplit(colnames(dds.rtx), '-'), function(r) paste(r[1], r[2], r[3], sep='-'))
stopifnot(all(names(surv_clust) %in% clin.uvm$submitter_id))

clin.uvm$clust.HERV <- surv_clust[clin.uvm$submitter_id]

TCGAanalyze_survival(clin.uvm, legend="", filename="analysisUM/07-pca_survival.pdf",
                     "clust.HERV", xlim=c(0,2500), dpi=600,
                     color=c(cA,cB),
                     main = "TCGA Set\n UVM", width=3.5, height=5, conf.int=F, risk.table=F)

pdf('analysisUM/pca_small.pdf', height=3.5, width=3)
PCAtools::biplot(p.rtx,
                 colby = 'metastatic',
                 colkey = c("No"="#666666", "Yes"="#990000"),
                 pointSize = 1,
                 hline = 0,
                 vline = 0,
                 legendPosition = 'none',
                 lab=F,
                 subtitleLabSize=0, titleLabSize = 0,
                 borderWidth = 0.2,
                 axisLabSize = 6
)
dev.off()
