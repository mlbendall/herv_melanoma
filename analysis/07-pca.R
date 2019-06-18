#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)

library(PCAtools)

#--- Load DESeq objects
load(snakemake@input[["deseq_rdata"]])

pdf(snakemake@output[[1]], paper='letter')

#--- Transcriptome
p <- PCAtools::pca(assay(tform.tx), metadata=colData(tform.tx), removeVar = 0.10)

PCAtools::screeplot(p,
    components = getComponents(p, 1:30),
    title="SCREE plot (transcriptome)",
    hline=50, vline=12
) +
geom_text(aes(5, 50, label = '50% explained variation', vjust = -1))

PCAtools::biplot(p,
    colby = 'category',
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="PCA biplot (transcriptome)"
)

PCAtools::pairsplot(p, 
    components = getComponents(p, c(1:3)),
    colby = 'category',
    pointSize = 0.4,
    title='PCA pairs (transcriptome)'
)

PCAtools::plotloadings(p,
    rangeRetain = 0.01,
    labSize = 3.0,
    title = 'Loadings (transcriptome)',
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24,
    col = c('limegreen', 'black', 'red3'),
    shapeSizeRange = c(5,5),
    drawConnectors = TRUE
)


#--- Retrotranscriptome

p <- PCAtools::pca(assay(tform.rtx), metadata=colData(tform.rtx), removeVar = 0.50)

PCAtools::screeplot(p,
    components = getComponents(p, 1:30),
    title="SCREE plot (retrotranscriptome)",
    hline=50, vline=21
) +
geom_text(aes(5, 50, label = '50% explained variation', vjust = -1))

PCAtools::biplot(p,
    colby = 'category',
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="PCA biplot (retrotranscriptome)"
)

PCAtools::pairsplot(p, 
    components = getComponents(p, c(1:3)),
    colby = 'category',
    pointSize = 0.4,
    title='PCA pairs (retrotranscriptome)'
)

PCAtools::plotloadings(p,
    rangeRetain = 0.01,
    labSize = 3.0,
    title = 'Loadings (retrotranscriptome)',
    subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24,
    col = c('limegreen', 'black', 'red3'),
    shapeSizeRange = c(5,5),
    drawConnectors = TRUE
)

dev.off()