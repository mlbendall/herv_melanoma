#! /usr/bin/env Rscript
library(tidyverse)
library(DESeq2)
library(matrixStats)
library(glmnet)
library(Boruta)
library(UpSetR)
library(c060)


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
  load("analysis2/05-filter_counts.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  outdat <- 'analysis2/10-feature_selection.Rdata'
  out_pdf <- 'analysis2/10-feature_selection/plots.pdf'
  
}


################################################################################
### Store selected variables in this list
################################################################################
selected_vars <- list()

makePDF <- TRUE
if(makePDF) dir.create(dirname(out_pdf), showWarnings = F)
if(makePDF) pdf(out_pdf, paper='letter', onefile=T)

################################################################################
### Retrotranscriptome expression data
################################################################################
countDat <- filt.herv
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(samples)))
samples$clust.herv <- clust.df$clust.retro.k4

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, samples, ~1)
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

################################################################################
### LRT
################################################################################
p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat, samples, ~ clust.herv)
dds_lrt <- DESeq2::DESeq(dds_lrt, test="LRT", reduced = ~ 1, parallel=T)
res_lrt <- results(dds_lrt, alpha=p.cutoff)
sig_lrt <- subset(res_lrt, padj < p.cutoff & abs(log2FoldChange) > lfc.cutoff)

selected_vars$lrt <- rownames(sig_lrt)

################################################################################
### Variable filtering
################################################################################
ntop <- nrow(tform) # No filtering
vars <- rowVars(assay(tform))
mat <- assay(tform)[order(-vars)[1:ntop], ]

cat(sprintf('%d variables\n', nrow(mat)))

################################################################################
### Randomized LASSO regression with stabiliy selection
################################################################################
set.seed(12345)
res <- c060::stabpath(samples$clust.herv, t(mat), weakness=0.1, 
                      family="multinomial",type.multinomial="grouped")

plot(res)

res.sel <- c060::stabsel(res, error=0.05, pi_thr=0.6, type="pfer")
length(res.sel$stable)

selected_vars$lasso <- names(res.sel$stable)

################################################################################
### Boruta
################################################################################
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=samples$clust.herv, doTrace=2, ntree=1000, maxRuns=1000)
print(bor.orig)
bor.model <- TentativeRoughFix(bor.orig)
print(bor.model)

selected_vars$boruta <- names(bor.model$finalDecision)[bor.model$finalDecision == 'Confirmed']


################################################################################
### Upset
################################################################################
# upset(fromList(selected_vars), sets=c('lasso', 'lassoU','boruta'), order.by = "freq")
upset(fromList(selected_vars), order.by = "freq")

# get the intersections
binmat <- fromList(selected_vars)
rn <- do.call(c, selected_vars)
rn <- rn[!duplicated(rn)]
rownames(binmat) <- rn

################################################################################
### SVM
################################################################################
set.seed(12345)

for(n in names(selected_vars)) {
  cat(sprintf('#--- %s - %d gene signature\n', n, sum(binmat[,n])))
  mat.sel <- mat[rownames(binmat[binmat[,n]==1,]), ]
  
  folds <- caret::createFolds(samples$clust.herv, k=5)
  svm.res <-lapply(folds, function(fld) {
    train <- t(mat.sel[,-fld])
    test <- t(mat.sel[,fld])
    train.lab <- samples[-fld, ]$clust.herv
    test.lab <- samples[fld, ]$clust.herv
    
    svm.model <- e1071::svm(y ~ . , data=data.frame(y=train.lab, train))
    pred <- predict(svm.model, data.frame(test))
    cm <- caret::confusionMatrix(table(pred, test.lab))
    list(svm.model, pred, cm)
  })
  
  svm.acc <- sapply(svm.res, function(l) l[[3]]$overall[1])
  cat(sprintf('Accuracy: %.5f\n', mean(svm.acc)))
}


dev.off()


binmat %>% rownames_to_column("rowid") %>% filter_at(vars(c("lrt", "lasso", "ss.lasso", "lassoU", "boruta")), ~.==1) %>% filter_at(vars(c("boruta")), ~.==0)
cat(sprintf('#--- %s - %d gene signature\n', n, sum(binmat[,n])))
mat.sel <- mat[rownames(binmat[binmat[,n]==1,]), ]

folds <- caret::createFolds(samples$clust.herv, k=5)
svm.res <-lapply(folds, function(fld) {
  train <- t(mat.sel[,-fld])
  test <- t(mat.sel[,fld])
  train.lab <- samples[-fld, ]$clust.herv
  test.lab <- samples[fld, ]$clust.herv
  
  svm.model <- e1071::svm(y ~ . , data=data.frame(y=train.lab, train))
  pred <- predict(svm.model, data.frame(test))
  cm <- caret::confusionMatrix(table(pred, test.lab))
  list(svm.model, pred, cm)
})

svm.acc <- sapply(svm.res, function(l) l[[3]]$overall[1])
cat(sprintf('Accuracy: %.5f\n', mean(svm.acc)))










mat.sel <- mat[rownames(binmat[binmat$boruta==1,]), ]

folds <- caret::createFolds(samples$clust.herv, k=5)
svm.res <-lapply(folds, function(fld) {
  train <- t(mat.sel[,-fld])
  test <- t(mat.sel[,fld])
  train.lab <- samples[-fld, ]$clust.herv
  test.lab <- samples[fld, ]$clust.herv
  
  svm.model <- e1071::svm(y ~ . , data=data.frame(y=train.lab, train))
  pred <- predict(svm.model, data.frame(test))
  cm <- caret::confusionMatrix(table(pred, test.lab))
  list(svm.model, pred, cm)
})

svm.acc <- sapply(svm.res, function(l) l[[3]]$overall[1])
mean(svm.acc)
















if(makePDF) dev.off()
save(selected_vars, partitions, 
     dds_lrt, 
     cv.fit.lasso, fit.lasso, 
     cv.fit.lassoU, fit.lassoU, 
     bor.model,
     binmat, inall, svm.res,
     file=outdat
     )



te.db[inall$rowid,] %>%
  select(Locus, Family, TE_CODING, TE_type, IntersectedGene) %>%
  write.table(sep='\t', quote=F)

