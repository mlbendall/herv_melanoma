library(tidyverse)
library(DESeq2)

load("analysis2/01-sample_data.Rdata")  
load("analysis2/02-clinical_data.Rdata")
load("analysis2/05-filter_counts.Rdata")
load("analysis2/31-additional_counts.Rdata")

out_dat <- 'analysis2/31a-study_pca.Rdata'

outdir <- 'analysis2/31-additional_counts'  
makePDF <- TRUE
if(makePDF) dir.create(outdir, showWarnings = F)


################################################################################
### Retrotranscriptome expression data
################################################################################
mfilt.herv.add <- counts.add.herv[rownames(mfilt.herv),]

countDat.comb <- cbind(mfilt.herv, mfilt.herv.add)

samples$study <- 'TCGA-UVM'
samples.comb <- rbind(samples, samples.add[,c(1:5,8)])

dds.comb <- DESeq2::DESeqDataSetFromMatrix(countDat.comb, samples.comb, ~1)
dds.comb <- DESeq2::DESeq(dds.comb, parallel=T)
tform.comb <- DESeq2::varianceStabilizingTransformation(dds.comb, blind=FALSE)

removeVar <- 0.3
pca.obj.comb <- PCAtools::pca(assay(tform.comb), metadata=samples.comb, removeVar=removeVar)
PCAtools::biplot(pca.obj.comb,
                 colby = 'study',
                 hline = 0,
                 vline = 0,
                 legendPosition = 'right',
                 lab=NULL
) + theme(aspect.ratio=1)

if(makePDF) ggsave(file.path(outdir, 'biplot.pdf'), paper='USr')

save(samples.comb, dds.comb, tform.comb, pca.obj.comb, file=out_dat)


################################################################################
### Batch adjustment
################################################################################
# adjusted <- sva::ComBat_seq(as.matrix(countDat.comb), batch=samples.comb$study)
# dds.adj <- DESeq2::DESeqDataSetFromMatrix(adjusted, samples.comb, ~1)
# dds.adj <- DESeq2::DESeq(dds.adj, parallel=T)
# tform.adj <- DESeq2::varianceStabilizingTransformation(dds.adj, blind=FALSE)
# 
# removeVar <- 0.3
# pca.obj.adj <- PCAtools::pca(assay(tform.adj), metadata=samples.comb, removeVar=removeVar)
# PCAtools::biplot(pca.obj.adj,
#                  colby = 'study',
#                  hline = 0,
#                  vline = 0,
#                  legendPosition = 'right',
#                  lab=NULL
# ) + theme(aspect.ratio=1)
# 
# if(makePDF) ggsave(file.path(outdir, 'biplot.pdf'), paper='USr')


# save(samples.comb, dds.comb, tform.comb, pca.obj.comb, adjusted, dds.adj, tform.adj, pca.obj.adj, file=out_dat)


