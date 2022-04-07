#! /usr/bin/env Rscript

library(tidyverse)
library(fpc)

#--- Load DESeq objects
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  load(snakemake@input[["samp_rdata"]])
  load(snakemake@input[["clin_rdata"]])
  load(snakemake@input[["clust_rdata"]])
  # out_dat <- snakemake@output[["rdata"]]
  outdir <- snakemake@output[["outdir"]]
} else {
  cat("Not using snakemake\n")
  load("analysis2/01-sample_data.Rdata")  
  load("analysis2/02-clinical_data.Rdata")
  load("analysis2/06-unsupervised.Rdata")
  outdat <- 'analysis2/08-cluststat.Rdata'
  outdir <- 'analysis2/08-cluststat'
}

makePDF <- TRUE
if(makePDF) dir.create(outdir, showWarnings = F)

cstats <- fpc::cluster.stats(
  d=dist(t(cDat), method = "euclidean"), 
  clustering=(as.numeric(clust.df$clust.retro.k4)+2)%%4 + 1,
  alt.clustering=as.numeric(mdata$clust.lncRNA)
)

# Rand index
clusterings <- names(mdata)[grep('^clust', names(mdata))]
rand.df <- lapply(names(clust.df), function(c1) {
  sapply(clusterings, function(c2) {
    cass1 <- as.numeric(clust.df[,c1])
    cass2 <- as.numeric(mdata[,c2])
    fossil::adj.rand.index(cass1, cass2)
  })
}) %>% bind_rows %>% data.frame
rownames(rand.df) <- names(clust.df)
  
rand.df %>%
  rownames_to_column("rclust") %>%
  tidyr::separate(rclust, c("x1","x2","k"), sep="\\.") %>%
  select(-c(x1,x2)) %>%
  tidyr::pivot_longer(2:7) %>%
  ggplot(aes(x=k, y=value, group=name, color=name)) + geom_line() + 
  labs(title="Cluster assignment similarity", ylab="Rand index") + theme_minimal()

if(makePDF) ggsave(file.path(outdir, 'randindex.pdf'), paper='USr')


  # fisher.test(table(mdata$gender, mdata$clust.SCNA))
  # fisher.test(table((ccp.obj[[4]]$consensusClass + 2) %% 4, ccp.obj[[4]]$consensusClass) )
  # fisher.test(table((ccp.obj[[4]]$consensusClass + 2) %% 4, ccp.obj[[4]]$consensusClass) )
  # 
  # fisher.test(table(clust.df$clust.retro.k4, mdata$clust.lncRNA), simulate.p.value = TRUE)
  # 
  # 
  # 
  # fisher.test(table(mdata$PLCB4, mdata$clust.methyl))
