#! /usr/bin/env Rscript

library(tidyverse)
library(xlsx)

################################################################################
### Conversion functions
################################################################################
# Convert numeric cluster to factor
clustfac <- function(x){
  lvl <- paste0('C', sort(unique(x)))
  return(factor(paste0('C', x), levels=lvl))
}

# Convert mutation data to factor
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


################################################################################
### File paths
################################################################################
if(exists("snakemake")) {
  cat("Running with snakemake\n")
  samp_rdata <- snakemake@input[["samp_rdata"]]
  table_xlsx <- snakemake@input[["tableS1"]]
  table_cib <- snakemake@input[["cibersort_table"]]
  out_dat <- snakemake@output[[1]]
  
} else {
  cat("Not using snakemake\n")
  samp_rdata <- "analysis2/01-sample_data.Rdata"
  table_xlsx <- "metadata/mmc2.xlsx"
  table_cib <- 'metadata/TCGA.Kallisto.fullIDs.cibersort.relative.tsv'
  out_dat <- "analysis2/02-clinical_data.Rdata"
}

load(samp_rdata)

################################################################################
### Load molecular-clinical data
################################################################################
molclin <- xlsx::read.xlsx(table_xlsx, 6, startRow=2, endRow=82, 
                           header=T, stringsAsFactors=F)
stopifnot(all(molclin$Patient.ID %in% samples$case_id))
stopifnot(all(samples$case_id %in% molclin$Patient.ID))
molclin$case_id <- factor(molclin$Patient.ID, levels=levels(samples$case_id))

samp.molclin <- dplyr::left_join(samples, molclin) %>%
  dplyr::mutate(
    death=factor(recode(Death..Metastasis, 
                        `Alive, no UM metastasis` = 'A',
                        `Alive, with UM metastasis` = 'AM',
                        `Death, metastatic UM` = 'DM',
                        `Death, other` = 'DO',
                        `Death, unknown` = 'DU'),
                 levels=c('A','AM','DM', 'DO', 'DU')),
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
  mutate(
    chr3CN = factor(recode(X3.CN..ABSOLUTE.,
                           `1`="1", 
                           `2, LOH`="2LOH",
                           `2, Subclonal LOH`="2SUB",
                           `2`="2",
                           `3`="3"),
                    levels=c("1","2LOH","2SUB","2","3")
    ),
    # Chr8CN = factor(X8q.CN..ABSOLUTE., levels=as.character(2:8)),
    chr8qCN = as.numeric(X8q.CN..ABSOLUTE.),
    chr8qISO = factor(Isochromosome.8q..ABSOLUTE., levels=c('Not Predicted','Predicted'))
  ) %>%
  mutate(
    D3M3 = factor(recode(chr3CN,
                         `1`="M3",
                         `2`="D3", 
                         .default="NA"),
                  levels=c("NA", "D3", "M3")
    )
  ) %>%
  select(-c(Death..Metastasis, SCNA.Cluster.No., DNA.Methyl.Cluster.No.,
            miRNA.Cluster.No., lncRNA.Cluster.No., mRNA.Cluster.No., 
            Paradigm.Cluster.No., X3.CN..ABSOLUTE., X8q.CN..ABSOLUTE.,
            Isochromosome.8q..ABSOLUTE.
  ))

# Mutated or not column
mutgenes <- c("GNAQ","GNA11","CYSLTR2","PLCB4","EIF1AX","SF3B1","SRSF2","BAP1")

tmp <- ifelse(samp.molclin[,mutgenes] == "NA", 0, 1)
colnames(tmp) <- paste0('mut.', colnames(tmp))

samp.molclin <- cbind(samp.molclin, tmp)
rm(tmp)


cat("------ Molecular/Clinical ------\n")
samp.molclin %>% head

################################################################################
### Load clinical-pathological data
################################################################################
clinpath <- xlsx::read.xlsx(table_xlsx, 5, startRow=2, endRow=82, 
                            header=T, stringsAsFactors=F)
stopifnot(all(clinpath$Patient.ID %in% samples$case_id))
stopifnot(all(samples$case_id %in% clinpath$Patient.ID))
clinpath$case_id <- factor(clinpath$Patient.ID, levels=levels(samples$case_id))


samp.clinpath <- dplyr::left_join(samples, clinpath) %>%
  dplyr::mutate(
    gender = factor(Gender, levels=c("Female", "Male")),
    stage = factor(gsub('^Stage ', '', AJCC.Clinical.Stage),
                   levels=c("IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")),
    tumor_diam_clin = suppressWarnings(as.numeric(Tumor.Basal.Diameter..Clinical.)),
    tumor_diam_path = suppressWarnings(as.numeric(Tumor.Basal.Diameter..pathological.)),
    # "Ciliary.Body..Pathology.Review."
    # "Ciliary.Body..TSS...Pathology.Data."
    # "Ciliary.Body..TSS.Data.Only."
    # "Histology.Cell.Type.Compare"
    tils_density=factor(TILS...Density, levels=c('Mild','Moderate','Heavy')),
    tams_density=factor(TAMS..Density, levels=c('Mild','Moderate','Heavy')),
    # "Closed.Connective.Loops"
    pigmentation=factor(Degree.of.Pigmentation, levels=c('Mild','Moderate','Heavy')),
    necrosis=factor(Necrosis, levels=c('Absent','Present')),
    # "Vital.Status"
    cause_of_death = factor(recode(Cause.of.Death,
                                   `[Not Applicable]`="NA", 
                                   `[Unknown]`="unknown",
                                   `Atrial fibrillation complications`='other',
                                   `Metastatic Pancreatic Cancer`='other',
                                   `Metastatic Uveal Melanoma`='MUM'
                                   ),
                            levels=c('NA','other','unknown','MUM')),
    metastatic = factor(recode(Metastatic.Disease,
                               Yes="Yes", No="No", YES="Yes"
                               ),
                        levels=c("No", "Yes")),
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
  select(-c(Gender, AJCC.Clinical.Stage, Tumor.Basal.Diameter..Clinical., 
            Tumor.Basal.Diameter..pathological., Degree.of.Pigmentation,
            Necrosis, Cause.of.Death, Metastatic.Disease,
            X1..Time.to.UM.Metastatsis, X2..Time.to.UM.Met.or.Death, 
            X3..Time.to.UM.Death.or.Last.Follow.up,
            TILS...Density, TAMS..Density
  ))

cat("------ Clinical/Pathology ------\n")
samp.clinpath %>% head

################################################################################
### Join clinical data
################################################################################
mdata <- left_join(samp.clinpath, samp.molclin)
row.names(mdata) <- mdata$sample_id

################################################################################
### Load leukocyte data
################################################################################
cat("----- Leukocyte fraction -------\n")
lf <- read.table('metadata/TCGA_all_leuk_estimate.masked.20170107.tsv',
                 sep='\t',
                 col.names=c('CancerType', 'SampleID', 'leukocyte_frac')
)
lf <- lf[lf$CancerType=='UVM',]
rownames(lf) <- sapply(strsplit(lf$SampleID, split="-"), function(x) paste(x[1], x[2], x[3], x[4], sep='-'))
lf <- lf[rownames(mdata),]
lf %>% head

#--- Join to mdata
stopifnot(all(rownames(mdata)==rownames(lf)))
mdata$leukocyte_frac <- lf$leukocyte_frac

################################################################################
### Load CIBERSORT data
################################################################################
cat("---------- CIBERSORT -----------\n")
cib <- read.table(table_cib, sep='\t', header=T)
cib <- cib[cib$CancerType=='UVM',]
rownames(cib) <- sapply(strsplit(cib$SampleID, split="\\."), function(x) paste(x[1], x[2], x[3], x[4], sep='-'))
stopifnot(all(rownames(mdata) %in% rownames(cib)))
cib <- cib[rownames(mdata), ]
cib %>% head

stopifnot(all(rownames(mdata)==rownames(cib)))
mdata <- cbind(mdata, cib[,3:24])

save(mdata, file=out_dat)
rm(list=setdiff(ls(), c('mdata')))

