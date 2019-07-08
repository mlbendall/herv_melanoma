#! /usr/bin/env Rscript

library(tidyverse)
library(fgsea)

# load('analysis/05-deseq.Rdata')
#--- Load DESeq objects
load(snakemake@input[["deseq_rdata"]])

#--- Load gene data
load(snakemake@input[["gene_rdata"]])

cat("--------- Gene to Sym ----------\n")
gsym %>% head

#--- Create ranked results

# adds symbol, contains duplicates
withdup <- lapply(results.tx, function(res) {
    res %>% 
        data.frame %>%
        tibble::rownames_to_column('row') %>%
        dplyr::inner_join(gsym, by=c("row"="GENEID")) %>% 
        as_tibble
})

# take mean of duplicates
rank.tx <- lapply(withdup, function(res){
    res %>%
        dplyr::select(SYM, stat) %>%
        na.omit() %>%
        distinct %>%
        dplyr::group_by(SYM) %>% 
        summarize(stat=mean(stat))
})
names(rank.tx) <- names(results.tx)

pathways <- list(
    hallmark = gmtPathways(snakemake@input[['h_gmt']]),
    oncogenic = gmtPathways(snakemake@input[['c6_gmt']]),
    immunologic = gmtPathways(snakemake@input[['c7_gmt']])  
)





  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
  

# Load the pathways into a named list
pathways.hallmark <- gmtPathways(snakemake@input[["deseq_rdata"]])




ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


topPathwaysUp <- fgseaRes[ES > 0, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)





pdf(snakemake@output[["pdf"]], paper='letter')


dev.off()

save(pathways, rank.tx, res_sym, results.fgsea, file=snakemake@output[["rdata"]])