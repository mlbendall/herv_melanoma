library(ggsci)

na_col <- "#dbdbda"
pri_mut <- '#b2586f'
sec_mut <- '#639c9c'
other_mut <- '#00bff6'

group_cols <- c(ggsci::pal_npg("nrc")(10), ggsci::pal_d3()(8))
names(group_cols) <- c("HERVK","HERVH","HERVW","HERVE","HERVF","HERVI","HERVL","ERVL","ERV1","HERVP","HERVS","PAB","MER4","PRIMA","ERV3","HARL","HUERS","MISC")

my_pals <- list(
  clust.mRNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
  clust.lncRNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
  clust.SCNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
  clust.miRNA = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6', 'C5'='#00ca77'),
  clust.paradigm = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6', 'C5'='#00ca77'),
  clust.methyl = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#ba3400', 'C4'='#0091c6'),
  clust.retro.k4 = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6'),
  clust.retro.k5 = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6', 'C5'='#00ca77'),
  clust.retro.k6 = c('C1'='#ffb769', 'C2'='#ff0000', 'C3'='#81d1e6', 'C4'='#0091c6',  'C5'='#00ca77', 'C6'="#df0192"),
  npg = ggsci::pal_npg("nrc")(10),
  stage = c("IIA"="#ebebeb", "IIB"="#dbdbdb", "IIIA"="#ababab", "IIIB"="#9b9b9b", "IIIC"="#8b8b8b", "IV"="#5b5b5b"),
  gender = c("Female"=na_col, "Male"="#000000"),
  death = c('A'='#f2f2f2', 'AM'='#919191', 'DM'='#000000', 'DO'='#bfbfbf', 'DU'='#bfbfbf'),
  metastatic = c("No"=na_col, "Yes"="#000000"),
  chr3CN=c("1"="#5ab9db", "2LOH"="#00a5ce", "2SUB"="#0087c1", "2"="#ffeee3", "3"="#ffd6ac"),
  chr8qCN=c("2"="#ffebdc", "3"="#ffcd9c", "4"="#ffa961", "5"="#ff8524", "6"="#ff5c00", "7"="#fd3800", "8"="#a52200"),
  chr8qISO=c("Not predicted"=na_col, "Predicted"="#a9a9a9"),
  BAP1 = c("NA"=na_col, "NON"="#a60e8b", "MIS"="#ff77ad", "FS"="#df0192", "SPL"="#8981c0", "IDEL"="#ffbec3", "HDEL"="#ffefff"),
  GNA11 = c("NA"=na_col, "R166H"=other_mut, "R183C"=other_mut, "Q209L"=pri_mut),
  GNAQ = c("NA"=na_col, "G48L"=sec_mut, "R183Q"=other_mut, 'Q209P'=pri_mut, 'Q209L'=pri_mut),
  CYSLTR2 = c("NA"=na_col, 'L129Q'=pri_mut),
  PLCB4 = c("NA"=na_col, 'D630N'=pri_mut),
  SF3B1 = c("NA"=na_col, "R625H"=pri_mut, "R625C"=pri_mut, "H662R"=other_mut, "T663P"=other_mut, "K666T"=sec_mut),
  SRSF2 = c("NA"=na_col, "Y92_H99del"=pri_mut, "Y92_H100del"=pri_mut, "S174_S179del"=pri_mut),
  EIF1AX = c("NA"=na_col, "G6D"=pri_mut, "K7_G8delinsR"=pri_mut, "G8R"=pri_mut, "G9R"=pri_mut, "G9D"=pri_mut, "G15D"=pri_mut, "W70R"=pri_mut),
  D3M3 = c("NA"=na_col, 'D3'="#ff0000", 'M3'="#0091c6"),
  herv.group = group_cols
)
names(my_pals[['npg']]) <- sprintf('category%d', 1:10)

rm(na_col, pri_mut, sec_mut, other_mut, group_cols)


plotCountsLine <- function (dds, gene, xaxis, facet=NULL, color=NULL, dfun=NULL, cisize=0.95,
                            normalized = TRUE, transform = TRUE, replaced = FALSE,
                            linecol = "red") {
  
  d <- DESeq2::plotCounts(dds, gene=gene, intgroup=names(colData(dds)), 
                          normalized=normalized, transform=transform, replaced=replaced,
                          returnData=TRUE)
  if(!is.null(dfun)) d <- d %>% dfun
  
  # Summary statistics
  d.ss <- d %>% group_by(.dots=lapply(c(xaxis, facet), as.symbol)) %>%
    summarise(N=length(count), avg=mean(count), sd=sd(count)) %>%
    mutate(se = sd / sqrt(N) ) %>%
    mutate(ci = se * qt(cisize/2 + 0.5, N-1))
  
  p <- ggplot(d, aes_string(x=xaxis, y="count"), inherit.aes = F)
  
  if(!is.null(color)) {
    p <- p + geom_point(aes_string(col=color), position=position_jitter(w=0.2,h=0), alpha=0.6)
  } else {
    p <- p + geom_point(position=position_jitter(w=0.2,h=0), alpha=0.6)
  }
  p <- p + 
    scale_y_log10(breaks=c(1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000),
                  minor_breaks=c(seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,5000,1000),seq(1000,5000,1000))) +
    labs(title=paste0(gene)) +
    geom_point(data=d.ss, inherit.aes = F,
               aes_string(x=xaxis, y="avg", group=1),
               colour=linecol) +
    geom_line(data=d.ss, inherit.aes = F,
              aes_string(x=xaxis, y="avg", group=1),
              colour=linecol)
  if(!is.null(facet)) p <- p + facet_grid(reformulate(facet))
  p
}


plotCountsErrorbar <- function (dds, gene, xaxis, facet=NULL, color=NULL, dfun=NULL, cisize=0.95,
                                normalized = TRUE, transform = TRUE, replaced = FALSE) {
  
  d <- DESeq2::plotCounts(dds, gene=gene, intgroup=names(colData(dds)), 
                          normalized=normalized, transform=transform, replaced=replaced,
                          returnData=TRUE)
  if(!is.null(dfun)) d <- d %>% dfun
  
  # Summary statistics
  d.ss <- d %>% group_by(.dots=lapply(c(xaxis, facet), as.symbol)) %>%
    summarise(N=length(count), avg=mean(count), sd=sd(count)) %>%
    mutate(se = sd / sqrt(N) ) %>%
    mutate(ci = se * qt(cisize/2 + 0.5, N-1))
  
  p <- ggplot(d, aes_string(x=xaxis, y="count"), inherit.aes = F)
  
  if(!is.null(color)) {
    p <- p + geom_point(aes_string(col=color), position=position_jitter(w=0.2,h=0), alpha=0.6)
  } else {
    p <- p + geom_point(position=position_jitter(w=0.2,h=0), alpha=0.6)
  }
  p <- p + 
    scale_y_log10(breaks=c(1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000),
                  minor_breaks=c(seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,5000,1000),seq(1000,5000,1000))) +
    labs(title=paste0(gene)) +
    geom_errorbar(data=d.ss, inherit.aes = F,
                  aes_string(x=xaxis, ymin="avg", y="avg", ymax="avg", group=xaxis),
                  width = 0.3, colour="red", alpha=0.8) +
    geom_errorbar(data=d.ss, inherit.aes = F,
                  aes_string(x=xaxis, ymin="avg-se", y="avg", ymax="avg+se", group=xaxis),
                  width = 0.3, colour="black", alpha=0.4)
  
  if(!is.null(facet)) p <- p + facet_grid(reformulate(facet))
  p
}


if(FALSE) {
  pdf('analysis2/legends.pdf', paper='USr')
  for(n in names(my_pals)) {
    plot.new()
    legend(0, 1, legend = names(my_pals[[n]]), 
           pch = 15, pt.cex = 3, cex = 1.5, bty = 'n',
           inset = c(-0.1, 0),
           title = n, 
           col = my_pals[[n]])
  }
  dev.off()
}



