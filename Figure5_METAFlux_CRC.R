require(dplyr)
require(tibble)
require(tidyr)
require(dplyr)
require(stringr)
require(readr)
require(readxl)
require(writexl)
require(ggplot2)
require(ggpubr)

# install.packages('./R_pkgs/spatstat.core_2.4-4.tar.gz', repos=NULL, type='source')
# devtools::install_github('KChen-lab/METAFlux')
library(METAFlux)
# for clustring
# install.packages('NbClust')

bulk_dat <- read_rds('~/Desktop/CRC_Macrophage_Figure5_6/TCGA-COAD.tpm_clinical.rds') %>% 
  as.matrix() %>% 
  t() #行为基因列为样本的矩阵
dim(bulk_dat)
# 肿瘤样本
bulk_dat <- bulk_dat[, colnames(bulk_dat)[str_detect(colnames(bulk_dat), '-01')]]
dim(bulk_dat)
# [1] 20052   469
data("human_gem")
data("human_blood")# for human derived samples
# calculate a single sample normalized MRAS from gene expression using GPR(Gene-protein-reaction).
scores<-calculate_reaction_score(bulk_dat)
# calculate the metabolic fluxes for the 13082 reactions.
flux<-compute_flux(mras=scores,medium = human_blood)
dim(flux)
# [1] 13082   469
write_rds(flux, file.path('~/Desktop/CRC_Macrophage_Figure5_6/METAFlux_flux_469.rds'))

#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux)){
    activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))
dim(all_pathway_score)
# 142 469
colnames(all_pathway_score) <- colnames(flux)
write_rds(all_pathway_score, '~/Desktop/CRC_Macrophage_Figure5_6/METAFlux_pathway_142.rds')
