library(tidyverse)
library(rpart)
library(rpart.plot)
library(DESeq2)
library(cowplot)

load("analysis2/10-feature_selection.Rdata")
source('analysis2/constants.R')
set.seed(12345)

################################################################################
### Get SVM from original TCGA data 
################################################################################
mat.sel <- t(mat[selected_vars$lasso,])
stopifnot(all(rownames(mat.sel) == rownames(samples)))
df.labels <- data.frame(y=samples$clust.herv, mat.sel)

rpart.fit <- rpart(y ~ . , data=df.labels, method = 'class')
rpart.plot(rpart.fit, type=5, extra=1, digits=3, box.palette = as.list(my_pals$clust.retro.k4))

svm.model <- e1071::svm(y ~ . , data=df.labels, probability=TRUE)

################################################################################
### Prediction
################################################################################
load("analysis2/31a-study_pca.Rdata")

# add cluster data to combined samples
is_tcga <- samples.comb$study == 'TCGA-UVM'
stopifnot(all(rownames(samples.comb[is_tcga,]) == rownames(samples)))
samples.comb$clust.herv <- NA
samples.comb[is_tcga,]$clust.herv <- as.character(samples$clust.herv)
samples.comb$clust.herv <- factor(samples.comb$clust.herv, levels=c('C1','C2','C3','C4'))

sel.comb <- t(assay(tform.comb)[selected_vars$lasso,])
pred <- predict(svm.model, data.frame(sel.comb))

# predicted cluster
samples.comb$pred.clust <- predict(svm.model, data.frame(sel.comb))
# cluster probabilities
tmp <- data.frame(attr(predict(svm.model, data.frame(sel.comb), probability=TRUE), 'probabilities'))
samples.comb$pred.C1 <- tmp$C1
samples.comb$pred.C2 <- tmp$C2
samples.comb$pred.C3 <- tmp$C3
samples.comb$pred.C4 <- tmp$C4

s.orig <- samples.comb[is_tcga,]
s.add <- samples.comb[!is_tcga,]

top_probs.orig <- c(s.orig[s.orig$pred.clust=='C1', ]$pred.C1,
               s.orig[s.orig$pred.clust=='C2', ]$pred.C2,
               s.orig[s.orig$pred.clust=='C3', ]$pred.C3,
               s.orig[s.orig$pred.clust=='C4', ]$pred.C4)

top_probs.add <- c(s.add[s.add$pred.clust=='C1', ]$pred.C1,
               s.add[s.add$pred.clust=='C2', ]$pred.C2,
               s.add[s.add$pred.clust=='C3', ]$pred.C3,
               s.add[s.add$pred.clust=='C4', ]$pred.C4)

# s.add %>% 
#   rownames_to_column('id') %>%
#   mutate(id=factor(id, levels=rev(id))) %>%
#   select(id, study, pred.C1, pred.C2, pred.C3, pred.C4) %>%
#   tidyr::gather(Clust, Prob, -c(id, study)) %>%
#   ggplot(aes(x=id, y=Prob, fill=Clust, group=study)) + 
#   geom_bar(stat="identity") + 
#   coord_flip()

plotCountsErrorbar(dds.comb, 'LTR46_Xq11.1', xaxis='study')

table(s.add[s.add$case_id=='Mel202',]$pred.clust)
table(s.add[s.add$case_id=='Mel270',]$pred.clust)
table(s.add[s.add$case_id=='92.1',]$pred.clust)
table(s.add[s.add$case_id=='OMM2.3',]$pred.clust)
table(s.add[s.add$case_id=='OMM1.3',]$pred.clust)

s.add[s.add$study=='GSE124720',]
s.add[s.add$study=='GSE149920',]
s.add[s.add$study=='GSE152705',]
s.add[s.add$study=='GSE176345',]
s.add[s.add$study=='GSE66719',]
s.add[s.add$study=='PRJNA596363',]

rpart.fit$variable.importance

pdf('analysis2/31-additional_counts/plotcounts.pdf', width=11, height=8, paper='USr')
grid.arrange(grobs=list(
plotCountsErrorbar(dds.comb, 'LTR46_Xq11.1', xaxis='study') + theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=1)),
plotCountsErrorbar(dds.comb, 'HML3_19q13.2', xaxis='study') + theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=1)),
plotCountsErrorbar(dds.comb, 'ERV316A3_6p21.33c', xaxis='study') + theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=1)),
plotCountsErrorbar(dds.comb, 'HERV3_19q13.42a', xaxis='study') + theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=1))
))
dev.off()



