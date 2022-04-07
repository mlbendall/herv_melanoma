#! /usr/bin/env Rscript

library(tidyverse)
library(xlsx)

# Convert numeric cluster to factor
clustfac <- function(x){
    lvl <- paste0('C', sort(unique(x)))
    return(factor(paste0('C', x), levels=lvl))
}


mutfac <- function(x){
    getlvl <- function(x) {
        mutpos <- function(x) {
            as.numeric(gsub("\\w(\\d+)\\w+", "\\1", x))
        }
        
        lvl <- unique(x)
        if("NA" %in% lvl) {
            lvl <- lvl[lvl != "NA"]
            lvl <- lvl[order(mutpos(lvl))]
            lvl <- c('NA', lvl)
        } else {
            lvl <- lvl[order(mutpos(lvl))]
            
        }
        lvl
    }
    return(factor(x, levels=getlvl(x)))
}

# Snakemake or not
if(exists("snakemake")) {
    cat("Running with snakemake\n")
    load(snakemake@input[["samp_rdata"]])
    table_xlsx <- snakemake@input[["tableS1"]]
    outfile <- snakemake@output[[1]]

} else {
    cat("Not using snakemake\n")    
    load("analysisUM/01-load_sample_data.Rdata")
    table_xlsx <- "metadata/mmc2.xlsx"
    outfile <- "analysisUM/01-load_clin_data.Rdata"
}

molclin <- xlsx::read.xlsx(table_xlsx, 6, startRow=2, endRow=82, 
                           header=T, stringsAsFactors=F)
stopifnot(all(molclin$Patient.ID %in% samples$case_id))
stopifnot(all(samples$case_id %in% molclin$Patient.ID))
molclin$case_id <- factor(molclin$Patient.ID, levels=levels(samples$case_id))

death_lvl = c('A','AM','DM', 'DO', 'DU')

samp.molclin <- dplyr::left_join(samples, molclin) %>%
    dplyr::mutate(
        death=factor(recode(Death..Metastasis, 
                            `Alive, no UM metastasis` = 'A',
                            `Alive, with UM metastasis` = 'AM',
                            `Death, other` = 'DO',
                            `Death, metastatic UM` = 'DM',
                            `Death, unknown` = 'DU'),
                     levels=death_lvl),
        clust.SCNA = clustfac(SCNA.Cluster.No.),
        clust.methyl = clustfac(DNA.Methyl.Cluster.No.),
        clust.miRNA = clustfac(miRNA.Cluster.No.),
        clust.lncRNA = clustfac(lncRNA.Cluster.No.),
        clust.mRNA = clustfac(mRNA.Cluster.No.),
        clust.paradigm = clustfac(Paradigm.Cluster.No.)
    ) %>%
    mutate(
        GNAQ = mutfac(GNAQ),
        GNA11 = mutfac(GNA11),
        CYSLTR2 = mutfac(CYSLTR2),
        PLCB4 = mutfac(PLCB4),
        EIF1AX = mutfac(EIF1AX),
        SF3B1 = mutfac(SF3B1),
        SRSF2 = mutfac(SRSF2),
        BAP1 = factor(recode(BAP1,
                             nonsense="NON",
                             missense_variant="MIS",
                             frameshift_variant="FS",
                             splicing_variant="SPL",
                             inframe_deletion="IDEL",
                             homozygous_deletion="HDEL"),
                      levels=c("NA","NON","MIS","FS","SPL","IDEL","HDEL")
        ),
    ) %>%
    select(-c(Death..Metastasis, SCNA.Cluster.No., DNA.Methyl.Cluster.No.,
              miRNA.Cluster.No., lncRNA.Cluster.No., mRNA.Cluster.No., 
              Paradigm.Cluster.No.
    ))

# Mutated or not column
mutgenes <- c("GNAQ","GNA11","CYSLTR2","PLCB4","EIF1AX","SF3B1","SRSF2","BAP1")

tmp <- ifelse(samp.molclin[,mutgenes] == "NA", 0, 1)
colnames(tmp) <- paste0('mut.', colnames(tmp))

samp.molclin <- cbind(samp.molclin, tmp)
rm(tmp)


cat("------ Molecular/Clinical ------\n")
samp.molclin %>% head


clinpath <- xlsx::read.xlsx(table_xlsx, 5, startRow=2, endRow=82, 
                            header=T, stringsAsFactors=F)
stopifnot(all(clinpath$Patient.ID %in% samples$case_id))
stopifnot(all(samples$case_id %in% clinpath$Patient.ID))
clinpath$case_id <- factor(clinpath$Patient.ID, levels=levels(samples$case_id))

yn_lvl = c("No", "Yes")
gender_lvl = c("Female", "Male")
stage_lvl = c("IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")

samp.clinpath <- dplyr::left_join(samples, clinpath) %>%
    dplyr::mutate(
        gender = factor(Gender, levels=gender_lvl),
        stage = factor(gsub('^Stage ', '', AJCC.Clinical.Stage),
                       levels=stage_lvl),
        metastatic = factor(recode(Metastatic.Disease,
                                   Yes="Yes",
                                   No="No",
                                   YES="Yes"
                                   ),
                            levels=yn_lvl),
        time_to_met = suppressWarnings(as.numeric(X1..Time.to.UM.Metastatsis)),
        time_to_met_or_death = suppressWarnings(as.numeric(X2..Time.to.UM.Met.or.Death)),
        time_to_death = suppressWarnings(as.numeric(X3..Time.to.UM.Death.or.Last.Follow.up))
    ) %>%
    mutate(
        stage_sim = factor(recode(stage,
                                  IIA="II", IIB="II", IIIA="III", IIIB="III", 
                                  IIIC="III"),
                           levels=c("II","III","IV"))
    ) %>%
    select(-c(Gender, AJCC.Clinical.Stage, Metastatic.Disease, 
              X1..Time.to.UM.Metastatsis, X2..Time.to.UM.Met.or.Death, 
              X3..Time.to.UM.Death.or.Last.Follow.up
    ))


cat("------ Clinical/Pathology ------\n")
samp.clinpath %>% head

mdata <- left_join(samp.clinpath, samp.molclin)
row.names(mdata) <- mdata$sample_id


cat("---------- CIBERSORT -----------\n")
cib <- read.table('metadata/TCGA.Kallisto.fullIDs.cibersort.relative.tsv', sep='\t', header=T)
cib <- cib[cib$CancerType=='UVM',]
rownames(cib) <- sapply(strsplit(cib$SampleID, split="\\."), function(x) paste(x[1], x[2], x[3], x[4], sep='-'))
stopifnot(all(rownames(mdata) %in% rownames(cib)))
cib <- cib[rownames(mdata), ]
cib %>% head

save(mdata, cib, file=outfile)
