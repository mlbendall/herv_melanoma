#! /usr/bin/env Rscript

library(tidyverse)
library(DESeq2)

library(PCAtools)

#--- Load DESeq objects
load(snakemake@input[["deseq_rdata"]])
load(snakemake@input[["clin_rdata"]])

stopifnot(all(colData(tform.rtx)$sample_id == colData(tform.tx)$sample_id))
stopifnot(all(colData(tform.rtx)$sample_id == row.names(mdata)))

pdf(snakemake@output[[1]], paper='letter')

#--- Transcriptome
p.tx <- PCAtools::pca(assay(tform.tx), metadata=mdata, removeVar = 0.10)
p.rtx <- PCAtools::pca(assay(tform.rtx), metadata=mdata, removeVar = 0.50)

color.by <- 'clust.SCNA'
colkey <- c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6')

PCAtools::screeplot(p.tx,
    components = getComponents(p.tx, 1:30),
    title="SCREE plot (transcriptome)",
    hline=50, vline=9
) +
geom_text(aes(5, 50, label = '50% explained variation', vjust = -1))

PCAtools::biplot(p.tx,
    colby = 'clust.SCNA',
    colkey = colkey,
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="PCA biplot (transcriptome)"
)

PCAtools::pairsplot(p.tx, 
    components = getComponents(p.tx, c(1:3)),
    colby = color.by,
    colkey = colkey,    
    pointSize = 0.4,
    title='PCA pairs (transcriptome)'
)

PCAtools::plotloadings(p.tx,
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
PCAtools::screeplot(p.rtx,
    components = getComponents(p.rtx, 1:30),
    title="SCREE plot (retrotranscriptome)",
    hline=50, vline=16
) +
geom_text(aes(5, 50, label = '50% explained variation', vjust = -1))

PCAtools::biplot(p.rtx,
    colby = color.by,
    colkey = colkey,    
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="PCA biplot (retrotranscriptome)"
)

PCAtools::pairsplot(p.rtx, 
    components = getComponents(p.rtx, c(1:3)),
    colby = color.by,
    colkey = colkey,
    pointSize = 0.4,
    title='PCA pairs (retrotranscriptome)'
)

PCAtools::plotloadings(p.rtx,
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

# Some other colorings
PCAtools::biplot(p.rtx,
    colby = "metastatic",
    colkey = c("No"="#666666", "Yes"="#990000"),
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="Retrotranscriptome, metastatic"
)

PCAtools::biplot(p.rtx,
    colby = "stage",
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="Retrotranscriptome, Stage"
)

PCAtools::biplot(p.rtx,
    colby = "stage_sim",
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="Retrotranscriptome, Stage (simplified)"
)

PCAtools::biplot(p.rtx,
    colby = "BAP1",
    colkey = c("NA"="#CCCCCC", "NON"="#a60e8b", "MIS"="#ff77ad", "FS"="#df0192", "SPL"="#8981c0", "IDEL"="#ffbec3", "HDEL"="#ffefff"),
    hline = 0,
    vline = 0,
    legendPosition = 'right',
    lab=F,
    title="Retrotranscriptome, BAP1"
)

dev.off()