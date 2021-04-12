" R script for Estrela et al (2021)
Nutrient dominance governs the assembly of microbial communities in mixed nutrient environments.

1) calculate Number of ESVs for each single and mixed CS treatment (S,A,SA,SS,AA)
2) calculate % of each family for each treatment
"

#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(reshape2)
library(tidyr)

#.. Assign path to data
#. usern = path to user directory
path.d<-paste0(usern, "/Estrelaetal_2021/data")
#.. Assign path to plots output
path.p<- paste0(usern,"/Estrelaetal_2021/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/Estrelaetal_2021/data/data_figures/")

#.. Read data 
setwd(path.d)
dat<-read.csv('data_master/data_16s_OD_singleCS_duoCS_processed.csv')

#.. 1) ---------------------------------------------------------------------------
#. determine Richness (Family and ESV numbers) for each community
min.ab<-0.0
dfp<-subset(dat, Relative_Abundance>min.ab)
ntaxa<- dfp %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS) %>%
  dplyr::summarise(
    nfam = length(unique(Family)),
    ngen = length(unique(Genus)),
    nesv = length(unique(ESV)))
ntaxa<-as.data.frame(ntaxa)

#. calculate number of ESVs (min and max) in each replicate communities for cs_type
ntaxa.esv<- ntaxa %>%
  dplyr::group_by(s_perc,cs_type,Experiment,nCS) %>%
  dplyr::summarise(
    nesv.min = min(nesv),
    nesv.max = max(nesv))
ntaxa.esv<-as.data.frame(ntaxa.esv)
ntaxa.esv

#.. 2) ---------------------------------------------------------------------------
#. sum abundance by Family for each community
min.ab<-0.00
dfp<-subset(dat, Relative_Abundance>min.ab)
sum_relab_fam<- dfp %>%
  dplyr::group_by(Inoculum,Carbon_Source,Replicate,s_perc,cs_type,Experiment,nCS,Family) %>%
  dplyr::summarise(
   sum_relab=sum(Relative_Abundance))
dfp<-as.data.frame(sum_relab_fam)

#. calculate % of each family for each CS type (S,A,SS,AA,SA) 
stat_summ<- dfp %>%
  dplyr::group_by(s_perc,cs_type,Experiment,nCS,Family) %>%
  dplyr::summarise(
    sum_relab_mean=round(mean(sum_relab),3),
    std = round(sd(sum_relab),3))
dfp.mean<-as.data.frame(stat_summ)

#. select dominant families
list.fam<-c("Enterobacteriaceae", "Pseudomonadaceae", "Moraxellaceae", "Rhizobiaceae")
dfp.s<-subset(select(dfp.mean,cs_type,Family,sum_relab_mean,std), Family %in% list.fam)
dfp.s

