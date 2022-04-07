#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(PCAtools)
library(ConsensusClusterPlus)

library(cowplot)
library(factoextra)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["clin_rdata"]])
  load(snakemake@input[["filt_rdata"]])
  out_dat <- snakemake@output[["rdata"]]
  outdir <- snakemake@output[["outdir"]]
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load("analysis2/05-filter_counts.Rdata")
  out_dat <- 'analysis2/06-unsupervised.Rdata'
  outdir <- 'analysis2/06-unsupervised'  
}

makePDF <- TRUE
if(makePDF) dir.create(outdir, showWarnings = F)

source('analysis2/constants.R')

################################################################################
### Retrotranscriptome expression data
################################################################################
countDat <- mfilt.herv
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(samples)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, samples, ~1)
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

################################################################################
### Principal component analysis
################################################################################
removeVar <- 0.3
pca.obj <- PCAtools::pca(assay(tform), metadata=mdata, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', removeVar*100, length(pca.obj$xvars)))

#--- SCREE plot
varline <- 50
varline.x <- min(which(cumsum(pca.obj$variance) >= varline))

horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
elbow <- PCAtools::findElbowPoint(pca.obj$variance)

PCAtools::screeplot(pca.obj,
                    components = getComponents(pca.obj, 1:30),
                    title="Retrotranscriptome SCREE",
                    hline=varline, vline=c(varline.x, horn$n, elbow)
) +
  geom_label(aes(x=varline.x+1, y=50, 
                 label = paste0(varline, '% var'), vjust = -1)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1))

if(makePDF) ggsave(file.path(outdir, 'scree_plot.pdf'), paper='USr', width=10, height=7.5)

cat(sprintf('%d PCs for Elbow method\n', elbow))
cat(sprintf('%d PCs for Horn method\n', horn$n))
cat(sprintf('%d PCs needed to explain %d percent of variation\n', varline.x, varline))

#--- Biplot (PC1 vs PC2)
colby <- 'clust.SCNA'
colkey <- my_pals[[colby]]
PCAtools::biplot(pca.obj,
                 colby = colby,
                 colkey = colkey,
                 shape='metastatic',
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("Retrotranscriptome PCA (", colby, ")")
) + theme(aspect.ratio=1)
if(makePDF) ggsave(file.path(outdir, 'biplot.pdf'), paper='USr')

#--- Pairs plot
maxcomp <- 5
PCAtools::pairsplot(pca.obj, 
                    components = getComponents(pca.obj, c(1:maxcomp)),
                    colby = colby,
                    colkey = colkey,
                    pointSize = 0.4,
                    title=paste0("Retrotranscriptome PC1-PC", maxcomp, " (", colby, ")")
) + theme(aspect.ratio=1)
if(makePDF) ggsave(file.path(outdir, 'pairs_plot.pdf'), paper='USr')


#--- Loadings
rangeRetain <- 0.01
PCAtools::plotloadings(pca.obj,
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
if(makePDF) ggsave(file.path(outdir, 'loadings.pdf'), paper='USr')


################################################################################
### Clustering
################################################################################
maxK <- 9
tDat <- assay(tform)[pca.obj$xvars, ]
cDat <- sweep(tDat, 1, apply(tDat, 1, median, na.rm=T))

ccp.obj <- ConsensusClusterPlus(cDat, maxK=maxK, reps=1000, pItem=0.8, pFeature=0.8,
                                title=outdir, clusterAlg="km", distance='euclidean',
                                seed=12345, plot="pdf")
icl <- calcICL(ccp.obj, title=outdir, plot='pdf')

stopifnot(all(rownames(pca.obj$metadata) == names(ccp.obj[[2]]$consensusClass)))

clust.df <- lapply(2:maxK, function(k) {
  factor(paste0('C', ccp.obj[[k]]$consensusClass), levels=paste0('C', 1:k))
})%>% bind_cols() %>% data.frame
colnames(clust.df) <- paste0('clust.retro.k', 2:maxK)
rownames(clust.df) <- names(ccp.obj[[2]]$consensusClass)

pca.obj$metadata <- cbind(pca.obj$metadata, clust.df)

wrap <- function(nclu) {
  columnname <- paste0('clust.retro.k', nclu)
  colkey <- ccp.obj[[nclu]]$clrs[[3]]
  names(colkey) <- paste0('C', 1:nclu)
  
  gp <- PCAtools::biplot(pca.obj,
                         colby = columnname,
                         colkey = colkey,
                         ellipse = TRUE,
                         ellipseConf = 0.9,
                         shape='death',
                         hline = 0,
                         vline = 0,
                         legendPosition = 'right',
                         lab=NULL,
                         title=paste0("Retrotranscriptome PCA (", columnname, ")")
  ) + 
    xlim(-45, 45) + 
    ylim(-35, 35) + 
    scale_shape_manual(values=c("A"=1, "AM"=16, "DM"=15, "DO"=2, "DU"=5)) + 
    theme(aspect.ratio=1)
  print(gp)
  print(table(pca.obj$metadata[ ,columnname], pca.obj$metadata$metastatic))
}

if(makePDF) pdf(file.path(outdir, 'cluster_biplot.pdf'), paper='USr')
for(nclu in 2:maxK) {
  wrap(nclu)
}
if(makePDF) dev.off()


cplots <- lapply(2:maxK, function(nclu) {
  p1 <- factoextra::fviz_cluster(
    list(data=t(cDat), cluster=ccp.obj[[nclu]]$consensusClass),
    ellipse.type="norm", ellipse.level=0.9,
    geom="point",
    palette = "npg"
  ) + theme_cowplot()
  sil <- cluster::silhouette(ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean"))
  p2 <- factoextra::fviz_silhouette(
    cluster::silhouette(ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean")),
    palette = "npg"
  )
  return(list(p1,p2))
  })
cplots <- unlist(cplots, recursive = FALSE)


if(makePDF) pdf(file.path(outdir, 'cluster_analysis.pdf'), width=8.5, height=11)
gridExtra::marrangeGrob(grobs = cplots, layout_matrix=matrix(1:8,4,2,T))
if(makePDF) dev.off()

sil.score <- sapply(2:maxK, function(k) {
  mean(cluster::silhouette(ccp.obj[[k]]$consensusClass, dist(t(cDat), method = "euclidean"))[,3])
})
ggplot(data.frame(x=2:maxK, y=sil.score)) + geom_point(aes(x,y)) + geom_line(aes(x,y)) + 
  xlab("k (num clusters)") + ylab('average silhouette score') + theme_cowplot()
if(makePDF) ggsave(file.path(outdir, 'silhouette_score.pdf'), paper='USr')

# Identify observation with negative silhouette
# neg_sil_index <- which(sil[, "sil_width"] < 0)
# sil[neg_sil_index, , drop = FALSE]

#--- Eigencor plot
if(makePDF) pdf(file.path(outdir, 'eigencorplot.pdf'), paper='USr', width=10, height=7.5)
PCAtools::eigencorplot(pca.obj,
                       components = getComponents(pca.obj, 1:horn$n),
                       metavars = c('Age', 'gender', 'metastatic', 'death', 'stage', 'Total.muts', 'time_to_met',
                                    'leukocyte_frac',
                                    'clust.retro.k2','clust.retro.k3','clust.retro.k4','clust.retro.k5',
                                    'clust.SCNA', 'clust.mRNA', 'clust.lncRNA', 'clust.miRNA'),
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)

PCAtools::eigencorplot(pca.obj,
                       components = getComponents(pca.obj, 1:horn$n),
                       metavars = c("GNAQ","GNA11","CYSLTR2","PLCB4","EIF1AX","SF3B1","SRSF2","BAP1"),
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)

PCAtools::eigencorplot(pca.obj,
                       components = getComponents(pca.obj, 1:horn$n),
                       metavars = c("B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive",
                                    "T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper",
                                    "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting",
                                    "NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting",
                                    "Dendritic.cells.activated","Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils"),
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
if(makePDF) dev.off()


################################################################################
### Recoding clusters
################################################################################
### k=4
table(mdata$clust.mRNA, clust.df$clust.retro.k4)

clust.df <- clust.df %>% 
  mutate(
    clust.retro.k4 = factor(recode(clust.retro.k4, 
                                   `C1`="C2",
                                   `C2`="C1",
                                   `C3`="C4",
                                   `C4`="C3"),
                            levels=c('C1','C2','C3','C4')
    )
  )

table(mdata$clust.mRNA, clust.df$clust.retro.k4)

### k=5
table(clust.df$clust.retro.k4, clust.df$clust.retro.k5)
clust.df <- clust.df %>% 
  mutate(
    clust.retro.k5 = factor(recode(clust.retro.k5,
                                   `C1`="C2",
                                   `C2`="C1",
                                   `C3`="C4",
                                   `C4`="C3",
                                   `C5`="C5"),
                            levels=c('C1','C2','C3','C4', 'C5')
    )
  ) 
table(clust.df$clust.retro.k4, clust.df$clust.retro.k5)

### k=6
table(clust.df$clust.retro.k5, clust.df$clust.retro.k6)
clust.df <- clust.df %>% 
  mutate(
    clust.retro.k6 = factor(recode(clust.retro.k6, 
                                   `C1`="C6",
                                   `C2`="C1",
                                   `C3`="C4",
                                   `C4`="C3",
                                   `C5`="C2",
                                   `C6`="C5"),
                            levels=c('C1','C2','C3','C4', 'C5', 'C6')
    )
  )
table(clust.df$clust.retro.k5, clust.df$clust.retro.k6)
rownames(clust.df) <- rownames(mdata)

pca.obj$metadata$clust.retro.k4 <- clust.df$clust.retro.k4
pca.obj$metadata$clust.retro.k5 <- clust.df$clust.retro.k5
pca.obj$metadata$clust.retro.k6 <- clust.df$clust.retro.k6

clust.df$group.k4 <- factor(recode(clust.df$clust.retro.k4, `C1`="C1C2", `C2`="C1C2", `C3`="C3C4", `C4`="C3C4"),
                            levels=c("C1C2", "C3C4"))

if(makePDF) pdf(file.path(outdir, 'cluster_final.pdf'), paper='USr')

PCAtools::biplot(pca.obj,
                 colby = 'clust.retro.k4',
                 colkey = my_pals[['clust.retro.k4']],
                 ellipse = TRUE,
                 ellipseConf = 0.9,
                 shape='death',
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HERV clustering (k=4)")
) + 
  xlim(-45, 45) + 
  ylim(-30, 30) + 
  scale_shape_manual(values=c("A"=1, "AM"=16, "DM"=15, "DO"=2, "DU"=5)) + 
  theme(aspect.ratio=1)
  

PCAtools::biplot(pca.obj,
                 colby = 'clust.retro.k5',
                 colkey = my_pals[['clust.retro.k5']],
                 ellipse = TRUE,
                 ellipseConf = 0.9,
                 shape='death',
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HERV clustering (k=5)")
) + 
  xlim(-45, 45) + 
  ylim(-30, 30) + 
  scale_shape_manual(values=c("A"=1, "AM"=16, "DM"=15, "DO"=2, "DU"=5)) + 
  theme(aspect.ratio=1)

PCAtools::biplot(pca.obj,
                 colby = 'clust.retro.k6',
                 colkey = my_pals[['clust.retro.k6']],
                 ellipse = TRUE,
                 ellipseConf = 0.9,
                 shape='death',
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL,
                 title=paste0("HERV clustering (k=6)")
) + 
  xlim(-45, 45) + 
  ylim(-30, 30) + 
  scale_shape_manual(values=c("A"=1, "AM"=16, "DM"=15, "DO"=2, "DU"=5)) + 
  theme(aspect.ratio=1)

if(makePDF) dev.off()

print(table(pca.obj$metadata$clust.retro.k4, pca.obj$metadata$metastatic))
print(table(pca.obj$metadata$clust.retro.k4, pca.obj$metadata$death))

print(table(pca.obj$metadata$clust.retro.k5, pca.obj$metadata$metastatic))
print(table(pca.obj$metadata$clust.retro.k6, pca.obj$metadata$metastatic))


save(pca.obj, ccp.obj, tDat, cDat, clust.df, file=out_dat)
