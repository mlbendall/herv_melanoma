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
### Create data partitions
################################################################################
# set.seed(12345)
# partitions <- lapply(caret::createDataPartition(samples$clust.herv, p=0.7, times=100), function(part) {
#     ssplit <- 1:nrow(samples) %in% part
#     names(ssplit) <- rownames(samples)
#     ssplit
# })

################################################################################
### LASSO regression (grouped)
################################################################################
# lasso_vars <- list()
# for(i in 1:100) {
#   ssplit <- partitions[[i]]
#   train <- mat[,ssplit]
#   train.lab <- samples[ssplit,]$clust.herv
#   test <- mat[,!ssplit]
#   test.lab <- samples[!ssplit,]$clust.herv
#   
#   set.seed(12345)
#   
#   # select lambda
#   cv.fit.lasso <- cv.glmnet(t(train), train.lab, alpha = 1,
#                             family="multinomial",
#                             type.multinomial = "grouped"
#   )
#   plot(cv.fit.lasso)
#   
#   best_lam <- cv.fit.lasso$lambda.1se
#   
#   fit.lasso <- glmnet(t(train), train.lab, alpha = 1,
#                       family="multinomial", type.multinomial = "grouped"
#   )
#   
#   
#   plot(fit.lasso, label=T, xvar="lambda")
#   abline(v = log(best_lam), lty = 2)
#   abline(v = log(cv.fit.lasso$lambda.min), lty = 2)
#   
#   confusion.glmnet(fit.lasso, newx=t(test), newy=test.lab, s=best_lam)
#   assess.glmnet(fit.lasso, newx=t(test), newy=test.lab, s=best_lam)
#   
#   cf <- coef(fit.lasso, s=best_lam)
#   tmp <- lapply(cf, function(m) as.matrix(m)[which(m != 0),])
#   # Since we used type.multinomial = "grouped" nonzero coef are same in all levels
#   stopifnot(all(names(tmp[[1]]) == names(tmp[[2]])))
#   stopifnot(all(names(tmp[[1]]) == names(tmp[[3]])))
#   stopifnot(all(names(tmp[[1]]) == names(tmp[[4]])))
#   sel.vars <- names(tmp[[1]])
#   
#   sel.vars <- sel.vars[sel.vars != '(Intercept)']
#   lasso_vars[[ paste0('lasso', i) ]] <- sel.vars
#   rm(tmp, sel.vars)
# }
# # get the intersections
# l.binmat <- fromList(lasso_vars)
# rn <- do.call(c, lasso_vars)
# rn <- rn[!duplicated(rn)]
# rownames(l.binmat) <- rn
# 
# l.binmat <- l.binmat[order(rowSums(l.binmat), decreasing = T),]
# rownames(l.binmat)[order(rowSums(l.binmat), decreasing = T)]
# 
# selected_vars$lasso <- selected_vars$lasso1

##### stability selection
set.seed(12345)
res <- c060::stabpath(samples$clust.herv, t(mat), weakness=0.1, 
                      family="multinomial",type.multinomial="grouped")

plot(res)

res.sel <- c060::stabsel(res, error=0.05, pi_thr=0.6, type="pfer")
length(res.sel$stable)

selected_vars$ss.lasso <- names(res.sel$stable)





################################################################################
### LASSO regression (ungrouped)
################################################################################
# set.seed(12345)
# 
# # select lambda
# cv.fit.lassoU <- cv.glmnet(t(train), train.lab, alpha = 1,
#                     family="multinomial"
# )
# plot(cv.fit.lassoU)
# 
# best_lamU <- cv.fit.lassoU$lambda.1se
# 
# fit.lassoU <- glmnet(t(train), train.lab, alpha = 1,
#               family="multinomial"
# )
# 
# 
# plot(fit.lassoU, xvar="lambda")
# abline(v = log(best_lamU), lty = 2)
# abline(v = log(cv.fit.lassoU$lambda.min), lty = 2)
# 
# confusion.glmnet(fit.lassoU, newx=t(test), newy=test.lab, s=best_lamU)
# assess.glmnet(fit.lassoU, newx=t(test), newy=test.lab, s=best_lamU)
# 
# cfU <- coef(fit.lassoU, s=best_lamU)
# tmp <- lapply(cfU, function(m) as.matrix(m)[which(m != 0),])
# sel.varsU <- unique(unlist(lapply(tmp, names)))
# sel.varsU <- sel.varsU[sel.varsU != '(Intercept)']
# selected_vars$lassoU <- sel.varsU
# rm(tmp, sel.varsU)


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
upset(fromList(lasso_vars), order.by = "freq")
upset(fromList(selected_vars), sets=c('lrt', 'ss.lasso', 'boruta'), order.by = "freq")

# get the intersections
binmat <- fromList(selected_vars)
rn <- do.call(c, selected_vars)
rn <- rn[!duplicated(rn)]
rownames(binmat) <- rn

inall <- binmat %>% rownames_to_column("rowid") %>% filter_at(vars(c("boruta", "ss.lasso")), ~.==1)

stopifnot(all(sapply(inall$rowid, function(n){
  all(sapply(selected_vars, function(lst) 'MER4_1p21.3' %in% lst))
})))

# in1110 <- binmat %>% rownames_to_column("rowid") %>% filter_at(vars(c("lrt", "lasso", "lassoU")), ~.==1) %>% filter_at(vars(c("boruta")), ~.==0)
# sapply(in1110$rowid, function(n){
#   all(c(n %in% selected_vars[[1]],n %in% selected_vars[[2]],n %in% selected_vars[[3]],!n %in% selected_vars[[4]]))
# })
# 
# in1000 <- binmat %>% rownames_to_column("rowid") %>% filter_at(vars(c("lrt")), ~.==1) %>% filter_at(vars(c("lasso", "lassoU", "boruta")), ~.==0)
# all(sapply(in1000$rowid, function(n){
#   all(c(n %in% selected_vars[[1]], !n %in% selected_vars[[2]], !n %in% selected_vars[[3]], !n %in% selected_vars[[4]]))
# }))
# 
# nrow(binmat[rowSums(binmat) >= 3,])

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

