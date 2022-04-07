#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)
library(WGCNA)
options(stringsAsFactors=FALSE)
# library(PCAtools)
# library(RColorBrewer)
# library(pheatmap)
# library(TCGAbiolinks)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  ##
} else {
  cat("Not using snakemake\n")
  load("analysisUM/05-deseq_int.Rdata")
  load("analysisUM/01-load_clin_data.Rdata")
  load("analysisUM/01-load_gene_data.Rdata")  
}

# Combined
stopifnot(all(colnames(dds.rtx) == colnames(dds.tx)))
counts.comb <- rbind(counts(dds.tx, normalize=F, replace=F), counts(dds.rtx, normalize=F, replace=F))
dds <- DESeqDataSetFromMatrix(counts.comb, colData(dds.tx), ~1)
dds <- DESeq(dds, parallel=T)

# Omit low count rows
cutoff.count <- 10
cutoff.samp <- floor(ncol(dds) * 0.05)

cat(sprintf("Filtering: Keep genes with %d or more counts in %d or more samples.\n", cutoff.count, cutoff.samp))
cat(sprintf("%d genes retained (%d original)\n", sum(rowSums(counts(dds)>=cutoff.count)>=cutoff.samp), dim(dds)[1]))

dds <- dds[rowSums(counts(dds)>=cutoff.count)>=cutoff.samp, ]
datExpr0 <- assay(dds)

# Check for good genes
gsg <- goodSamplesGenes(datExpr0, verbose = 5)
stopifnot(gsg$allOK)

# sampleTree <- hclust(dist(t(datExpr0)), method='average')
# plot(sampleTree)

# Transform data
vsd <- getVarianceStabilizedData(dds)
datExpr <- t(vsd)

# WGCNA Analysis
maxPOutliers <- 0.05

powers <- c(c(1:10), seq(from = 12, to=40, by=2))
sft <- pickSoftThreshold(
  datExpr,
  corFnc = "bicor", 
  corOptions=list(maxPOutliers=maxPOutliers),
  networkType = "signed hybrid",
  powerVector = powers,
  verbose = 5
)


################################################################################
### Plot the scale independence
################################################################################
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <- 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=c(0.8, 0.85, 0.9), col="#99000099", lty=4)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


################################################################################
### The scale-free topology index exceeds 0.8 for soft power = 4,
### 0.85 for soft power = 12, and 0.9 for soft power = 22
### Using soft power = 7 since it is small peak
################################################################################

softpower <- 7

################################################################################
### Topological Overlap Matrix
################################################################################
# adjacency <- adjacency(datExpr, power=softpower)
# TOM <- TOMsimilarity(adjacency)
# dissTOM = 1 - TOM
# geneTree <- hclust(as.dist(dissTOM), method="average")
# 
# sizeGrWindow(12, 9)
# plot(geneTree, xlab="", sub="", main="Gene clustering based on TOM-based similarity",
#      labels=FALSE, hang=0.04)
# 
# minModuleSize <- 30
# dynamicMods <- cutreeDynamic(dendro = geneTree, 
#                              distM = dissTOM,
#                              deepSplit = 2,
#                              pamRespectsDendro = FALSE,
#                              minClusterSize = minModuleSize
#                              )
# table(dynamicMods)
# # There are 51 modules ranging from 7441 genes to 30 genes
# # There are 285 unassigned genes
# 
# dynamicColors <- labels2colors(dynamicMods)
# table(dynamicColors)
# sizeGrWindow(8, 6)
# plotDendroAndColors(geneTree,
#                     dynamicColors,
#                     "Dynamic Tree Cut",
#                     dendroLabels=FALSE,
#                     hang = 0.03,
#                     addGuide = TRUE,
#                     guideHang = 0.05,
#                     main = "Gene Dendrogram and module colors"
#                     )

################################################################################
### One-step
################################################################################
net <- blockwiseModules(
  # Input data
  datExpr, 
  # Options for splitting data into blocks
  maxBlockSize = 30000, ## 5000,
  # Network construction arguments: correlation options
  corType = "bicor", ## "pearson",
  maxPOutliers = maxPOutliers, ## 1, 
  # Adjacency function options
  power = softpower,  
  networkType = "signed hybrid", ## "unsigned",
  # Basic tree cut options
  minModuleSize = 30, ## min(20, ncol(datExpr)/2 ),
  # Advanced tree cut options
  pamRespectsDendro = FALSE, ## TRUE,
  # Options controlling behaviour
  numericLabels = TRUE,
  nThreads = 8,
  verbose = 3
)



table(net$colors)
table(net9$colors)


sizeGrWindow(12, 9)
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    main = "Gene Dendrogram and module colors",
                    dendroLabels=FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# softpower9 <- 22
# net9 <- blockwiseModules(
#   # Input data
#   datExpr, 
#   # Options for splitting data into blocks
#   maxBlockSize = 30000, ## 5000,
#   # Network construction arguments: correlation options
#   corType = "bicor", ## "pearson",
#   maxPOutliers = maxPOutliers, ## 1, 
#   # Adjacency function options
#   power = softpower9,  
#   networkType = "signed hybrid", ## "unsigned",
#   # Basic tree cut options
#   minModuleSize = 30, ## min(20, ncol(datExpr)/2 ),
#   # Advanced tree cut options
#   pamRespectsDendro = FALSE, ## TRUE,
#   # Options controlling behaviour
#   numericLabels = TRUE,
#   nThreads = 8,
#   verbose = 3
# )
# 
# table(net9$colors)
# 
# sizeGrWindow(12, 9)
# mergedColors9 <- labels2colors(net9$colors)
# plotDendroAndColors(net$dendrograms[[1]], 
#                     colors=data.frame(
#                       mergedColors[net$blockGenes[[1]]],
#                       mergedColors9[net9$blockGenes[[1]]]
#                     ),
#                     main = "Gene Dendrogram and module colors",
#                     dendroLabels=FALSE,
#                     hang = 0.03,
#                     addGuide = TRUE,
#                     guideHang = 0.05
# )



################################################################################
### Relating modules to clinical traits
################################################################################

nGenes <- ncol(datExpr)
nSamples <- ncol(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

datTraits <- data.frame(
  Age=mdata$Age,
  isFemale=ifelse(mdata$gender=='Female', 1, 0),
  isMetastatic=ifelse(mdata$metastatic=='Yes', 1, 0),
  
  clust.SCNA=as.numeric(mdata$clust.SCNA),
  clust.methyl=as.numeric(mdata$clust.methyl),
  clust.miRNA=as.numeric(mdata$clust.miRNA),
  clust.lncRNA=as.numeric(mdata$clust.lncRNA),
  clust.mRNA=as.numeric(mdata$clust.mRNA),
  clust.paradigm=as.numeric(mdata$clust.paradigm),
  BAP1=as.numeric(mdata$BAP1)
)

moduleTraitCor <- WGCNA::cor(MEs, as.matrix(datTraits), use="p")
moduleTraitPvalue <- WGCNA::corPvalueFisher(moduleTraitCor, nSamples)

sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix  <- paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

CorLevelPlot(datTraits, x=colnames(datTraits), y=colnames(datTraits))
CorLevelPlot(cbind(MEs, datTraits),  x=colnames(datTraits), y=colnames(MEs))
CorLevelPlot(cbind(MEs, datTraits),  x=colnames(datTraits), y=colnames(MEs))


CorLevelPlot(mdata,  x=colnames(mdata), y=colnames(mdata))

# cat("\n----------- Combined  ----------\n")
# summary(res)
# save(padj_cutoff, lfc_cutoff, dds, tform, res, file="analysis/05-deseq.Rdata")
