library(tidyverse)
library(fgsea)
library(DESeq2)

load('analysis2/03-gene_data.Rdata')
load('analysis2/09-de.Rdata')


if(FALSE) {
  tof <- lapply(retro.fam$family, function(f) {
    rownames(subset(retro.annot, family==f))
  })
  cat(retro.fam$family[1], "\tNA\t", paste(tof[[1]], collapse = '\t'), '\n', sep='', file="hervfam.gmt")
  for(i in 2:length(tof)){
    cat(retro.fam$family[i], "\tNA\t", paste(tof[[i]], collapse = '\t'), '\n', sep='', file="hervfam.gmt", append=T)
  }
  
  toc <- lapply(names(table(retro.annot$chrom))[1:24], function(f) {
    rownames(subset(retro.annot, class=='HERV' & chrom==f))
  })
  cat(names(table(retro.annot$chrom))[1], "\tNA\t", paste(toc[[1]], collapse = '\t'), '\n', sep='', file="hervchrom.gmt")
  for(i in 2:length(toc)){
    cat(names(table(retro.annot$chrom))[i], "\tNA\t", paste(toc[[i]], collapse = '\t'), '\n', sep='', file="hervchrom.gmt", append=T)
  }
}


#####

pathways <- fgsea::gmtPathways("hervchrom.gmt")

ct <- names(res)[5]
ranks <- res[[ct]][complete.cases(res[[ct]]),] %>% 
  data.frame %>%
  rownames_to_column('SYM') %>%
  select(SYM, stat) %>%
  deframe

res.fgsea <- fgseaMultilevel(pathways=pathways, stats=ranks, eps=1e-14) %>%
  as_tibble() %>%
  arrange(desc(NES))

res.fgsea %>%
  na.omit %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  scale_fill_manual(values=c("#CCCCCC", "#880000")) +
  coord_flip() +
  labs(x="chromosome", y="Normalized Enrichment Score",
       title=paste0('Chromosome - ', ct)) +
  theme_minimal()

ggsave('analysis2/09-de/chromGSEA.pdf', paper='letter')

#####


pathways <- fgsea::gmtPathways("hervfam.gmt")

ct <- names(res)[5]
ranks <- res[[ct]][complete.cases(res[[ct]]),] %>% 
  data.frame %>%
  rownames_to_column('SYM') %>%
  select(SYM, stat) %>%
  deframe

res.fgsea <- fgseaMultilevel(pathways=pathways, stats=ranks, eps=1e-14) %>%
  as_tibble() %>%
  arrange(desc(NES))


res.fgsea %>%
  na.omit %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  scale_fill_manual(values=c("#CCCCCC", "#880000")) +
  coord_flip() +
  labs(x="chromosome", y="Normalized Enrichment Score",
       title=paste0('HERV family - ', ct)) +
  theme_minimal()

ggsave('analysis2/09-de/familyGSEA.pdf', paper='letter')

pdf('analysis2/09-de/familyGSEA_enrich.pdf', width=7, height=3.5)
n <- 'HML2'
plotEnrichment(pathways[[n]], ranks) + labs(title=n)

n <- 'HERV9'
plotEnrichment(pathways[[n]], ranks) + labs(title=n)
dev.off()






plotEnrichment(pathways[["HERV9"]], ranks) + 
  labs(title=head(res.fgsea[order(res.fgsea$pval), ], 1)$pathway)
