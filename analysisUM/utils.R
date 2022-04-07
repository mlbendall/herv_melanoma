

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
